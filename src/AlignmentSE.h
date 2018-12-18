/*
 * AlignmentSE.h
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#ifndef SRC_ALIGNMENTSE_H_
#define SRC_ALIGNMENTSE_H_

#define _USE_MATH_DEFINES

#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <boost/algorithm/string/regex.hpp>
#include "DNAseq.h"
#include "ScoreScheme.h"
#include "MEM.h"
#include "MEMS.h"
#include "MetaGenome.h"
#include "BAM.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace Eigen;
using std::basic_string;
using std::istream;
using std::ostream;
using EGriceLab::SAMtools::BAM;
using std::vector;

/**
 * an alignment region of between a single-end read/query and a database/target
 */
struct AlignmentSE {
	typedef uint32_t CIGAR_OP_TYPE;
	/* nested types and enums */
	/** A SeedPair is a pair of exact matched seed region between query and target */
	struct SeedPair {
		/* constructors */
		/** default constructor */
		SeedPair() = default;

		/** construct from given values */
		SeedPair(uint32_t from, uint32_t start, uint32_t len, int32_t tid, uint32_t qShift = 0, uint32_t tShift = 0)
		: from(from), start(start), to(from + len), end(start + len), tid(tid), qShift(qShift), tShift(tShift)
		{ 	}

		/* member methods */
		uint32_t length() const {
			return end - start;
		}

		/* member fields */
		uint32_t from = 0; /* 0-based query position */
		uint32_t to = 0;   /* 1-based query position */
		uint32_t start = 0; /* 0-based target position */
		uint32_t end = 0; /* 1-based target position */
		uint32_t qShift = 0;  /* query shift */
		uint32_t tShift = 0;  /* target shift relative to chromosome */
		int32_t tid = -1; /* tid as chromId/locId determined by database */

		/**
		 * shift on target
		 * @param shift  positions to shift, positive indicating shifting right
		 */
		SeedPair& shiftTarget(int32_t shift) {
			start -= shift;
			end -= shift;
			tShift += shift;
			return *this;
		}

		/**
		 * shift on query
		 * @param shift  positions to shift, positive indicating shifting right
		 */
		SeedPair& shiftQuery(int32_t shift) {
			from -= shift;
			to -= shift;
			qShift += shift;
			return *this;
		}

		/** test whether two SeedPairs overlap (on target) */
		static bool isOverlap(const SeedPair& lhs, const SeedPair& rhs) {
			return lhs.tid == rhs.tid && lhs.start < rhs.end && lhs.end > rhs.start;
		}

		/** test whether two SeedPairs are compatitable (no overlaps) */
		static bool isCompatitable(const SeedPair& lhs, const SeedPair& rhs) {
			return lhs.tid == rhs.tid && !(lhs.start < rhs.end && lhs.end > rhs.start);
		}
	};

	/**
	 * a SeedMatch is an ordered lists of SeedPairs
	 */
	struct SeedMatch : public std::vector<SeedPair> {
		/* member methods */
		/** test whether this SeedMatch is compatitable */
		bool isCompatitable() const;

		/** get total length of SeedPairs within */
		uint32_t length() const;

		/** get tid of this SeedMatch */
		int32_t getTId() const {
			return front().tid;
		}

		/** get start position of this SeedMatch as the first SeedPair start */
		uint32_t getStart() const {
			return front().start;
		}

		/** get end position of this SeedMatch as the last SeedPair end */
		uint32_t getEnd() const {
			return back().end;
		}

		/** get start position of this SeedPair given max indel-rate */
		uint32_t getStart(double maxIndelRate) const {
			return getStart() - front().from * maxIndelRate;
		}

		/** get end position of this SeedPair given max indel-rate */
		uint32_t getEnd(uint32_t qLen, double maxIndelRate) const {
			return getEnd() + (qLen - back().to) * maxIndelRate;
		}

		/** shift SeedMatch on target */
		SeedMatch& shiftTarget(int32_t shift) {
			for(SeedMatch::value_type& seed : *this)
				seed.shiftTarget(shift);
			return *this;
		}

		/** shift SeedMatch on query */
		SeedMatch& shiftQuery(int32_t shift) {
			for(SeedMatch::value_type& seed : *this)
				seed.shiftQuery(shift);
			return *this;
		}

	};

	typedef vector<SeedMatch> SeedMatchList;

	typedef basic_string<uint32_t> state_str;
	/* constructors */
	/** default constructor */
	AlignmentSE() = default;

	/** construct an AlignmentSE with all fields */
	AlignmentSE(const DNAseq& query, const DNAseq& target,
			const string& qname, int32_t tid, const ScoreScheme* ss,
			const QualStr& qual, int32_t qShift, int32_t tShift,
			uint16_t flag, uint8_t mapQ, int32_t mtid, int32_t mpos, int32_t isize, uint32_t id)
	: query(query), target(target), qname(qname), tid(tid), ss(ss), qual(qual), qShift(qShift), tShift(tShift),
	  flag(flag), mapQ(mapQ), mtid(mtid), mpos(mpos), isize(isize), id(id),
	  qLen(query.length()), tLen(target.length()),
	  M(qLen + 1, tLen + 1), I(qLen + 1, tLen + 1), D(qLen + 1, tLen + 1) {
		if(this->qual.length() != query.length())
			this->qual.resize(qLen, QualStr::INVALID_Q_SCORE);
		initScores();
	}

	/** construct an AlignmentSE with required fields */
	AlignmentSE(const DNAseq& query, const DNAseq& target,
			const string& qname, int32_t tid, const ScoreScheme* ss,
			int32_t tShift = 0, const QualStr& qual = QualStr(), uint16_t flag = 0)
	: query(query), target(target), qname(qname), tid(tid), ss(ss), tShift(tShift), qual(qual), flag(flag),
	  qLen(query.length()), tLen(target.length()),
	  M(qLen + 1, tLen + 1), I(qLen + 1, tLen + 1), D(qLen + 1, tLen + 1) {
		if(this->qual.length() != query.length())
			this->qual.resize(qLen, QualStr::INVALID_Q_SCORE);
		initScores();
	}

	/* member methods */
	/** test whether this alignment is initiated */
	bool isInitiated() const {
		return !query.empty() && !target.empty() && !qual.empty() && ss != nullptr;
	}

	/** initiate all score matrices */
	AlignmentSE& initScores();

	/** calculate all scores in the entire region using Dynamic-Programming, return the alnScore as the maximum score found */
	AlignmentSE& calculateScores() {
		calculateScores(0, qLen, 0, tLen);
		alnScore = M.maxCoeff(&alnTo, &alnEnd); // determine aign 3' and score simultaneously
		return *this;
	}

	/** calculate all scores in the restricted "SeedMatch" regions using Dynamic-Programming, return the alnScore as the maximum score found */
	AlignmentSE& calculateScores(const SeedMatch& seeds);

	/** backtrace alnPath */
	AlignmentSE& backTrace();

	/** get 5' soft-clip size */
	uint32_t getClip5Len() const {
		return alnFrom - 0;
	}

	/** get 3' soft-clip size */
	uint32_t getClip3Len() const {
		return query.length() - alnTo;
	}

	/** get pos of */

	/** get align query length */
	uint32_t getAlnQLen() const {
		return alnTo - alnFrom;
	}

	/** get align target length */
	uint32_t getAlnTLen() const {
		return alnEnd - alnStart;
	}

	/** get align length, determined by alnPath, including M=XIDX but not H,P,N */
	uint32_t getAlnLen() const;

	/** get decoded aligned query seq */
	string getAlnQSeq() const;

	/** get decoded aligned target seq */
	string getAlnTSeq() const;

	/** get decoded aligned seq, alias to getAlnTSeq() */
	string getAlnSeq() const {
		return getAlnTSeq();
	}

	/** get cigar_str */
	BAM::cigar_str getAlnCigar() const;

	/** get MD:Z tag for this Alignment */
	string getAlnMDTag() const;

	/** get alignment score of this alignment, ignore soft-clips */
	double getAlnScore() const {
		return alnScore;
	}

	/** get overall score of this alignment, including soft-clips */
	double getScore() const {
		return alnScore - ss->clipPenalty * (getClip5Len() + getClip3Len());
	}

	/** evaluate this alignment log-liklihood using seq, align-path and quality */
	void evaluate();

	/** get log10-liklihood of this alignment */
	double log10lik() const {
		return log10P;
	}

	/** get log-liklihood of this alignment */
	double loglik() const {
		return log10lik() / ::log10(M_E);
	}

	/** get the p-value of this alignment */
	double pvalue() const {
		return ::pow(10.0, log10lik());
	}

	/**
	 * export core info of this alignment to a BAM record, with no aux data
	 */
	BAM exportBAM() const {
		return BAM(qname, flag, tid, tShift + alnStart, mapQ,
				getAlnCigar(), qLen, nt16Encode(query), qual,
				mtid, mpos, isize, id);
	}

	/* member fields */
	/* htslib required fields */
	DNAseq query;  // query
	DNAseq target; // target
	QualStr qual; // query qual
	string qname; // query name

	int32_t qShift = 0; // query shift relative to original read
	int32_t tShift = 0; // target shift relative to the target (chromo)
	int32_t tid = -1;  // target id, should be determined from the database
	const ScoreScheme* ss = nullptr;
	uint16_t flag = 0;
	uint8_t mapQ = INVALID_MAP_Q;
	int32_t mtid = -1;
	int32_t mpos = -1;
	int32_t isize = 0;
	uint64_t id = 0;

	uint32_t qLen = 0; // query length
	uint32_t tLen = 0; // target length

	MatrixXd M; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that x[i] is aligned to y[j]
	MatrixXd D; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that x[i] is aligned to a gap (deletion for query)
	MatrixXd I; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that y[j] is aligned to a gap (insertion for query)

	uint32_t alnFrom = 0; // 0-based alignment from
	uint32_t alnTo = 0; // 1-based alignment to
	uint32_t alnStart = 0; // 0-based alignment start
	uint32_t alnEnd = 0; // 1-based alignment end
	double alnScore = infV; // alignment score
	state_str alnPath; // backtrace alignment path using cigar ops defined htslib/sam.h
//	BAM::cigar_str alnCigar; // cigar_str compressed from alnPath
	double score = infV; // overall score of alignment, including 5'/3' clips
	double log10P = NAN;  // log10-liklihood of this alignment, given the alignment and quality
	double postP = NAN;   // posterior probability of this alignment, given all candidate alignments

	/* static fileds */
	static const uint8_t INVALID_MAP_Q = 0xff;
	static const uint32_t MAX_ITER = UINT16_MAX;
	static const double DEFAULT_INDEL_RATE;
	static const double MAX_INDEL_RATE;
	static const string STATES; // human readable STATES
	static const double DEFAULT_SCORE_REL_EPSILON; // default relative min score comparing to best score
	static const char ALIGN_GAP = '-';
	static const char ALIGN_MD_INS = '^'; // symbol use in MD:Z tag
	static const boost::regex MDTAG_LEADING_PATTERN;
	static const boost::regex MDTAG_MAIN_PATTERN;

	/* static methods */
	static CIGAR_OP_TYPE matchMax(double match, double ins, double del);
	static CIGAR_OP_TYPE insMax(double match, double ins);
	static CIGAR_OP_TYPE delMax(double match, double del);

	/* static methods */
	/**
	 * convert our DNAseq to BAM::seq_str encoding from htslib
	 */
	static BAM::seq_str nt16Encode(const DNAseq& seq);

	/** get algnQLen by MD tag */
	static uint32_t mdTag2alnQLen(const string& mdTag);

	/** get a search region given max indel-rate */
	static Loc getSearchRegion(uint64_t start, uint64_t end);

	/** get SeedMatchList from a pre-calculated MEMS and its shift of this region */
	static SeedMatchList getSeedMatchList(const MetaGenome& mtg, const MEMS& mems,
			uint32_t maxIt = MAX_ITER);

	/**
	 * calculate mapQ as posterior probability of candidate alignments using a uniform prior
	 * set the mapQ value of all alignments
	 */
	static vector<AlignmentSE>& calcMapQ(vector<AlignmentSE>& alnList);

	/* utility methods */
	/**
	 * calculate all scores in a given region using Dynamic-Programming
	 * @param from, start  0-based
	 * @param to, end  1-based
	 */
	void calculateScores(uint32_t from, uint32_t to, uint32_t start, uint32_t end);
};

inline AlignmentSE::CIGAR_OP_TYPE AlignmentSE::matchMax(double match, double ins, double del) {
	double max = match;
	int s = BAM_CMATCH;
	if(ins > max) {
		max = del;
		s = BAM_CINS;
	}
	if(del > max) {
//		min = del;
		s = BAM_CDEL;
	}
	return s;
}

inline AlignmentSE::CIGAR_OP_TYPE AlignmentSE::insMax(double match, double ins) {
	return match > ins ? BAM_CMATCH : BAM_CINS;
}

inline AlignmentSE::CIGAR_OP_TYPE AlignmentSE::delMax(double match, double del) {
	return match > del ? BAM_CMATCH : BAM_CDEL;
}

inline void AlignmentSE::calculateScores(uint32_t from, uint32_t to, uint32_t start, uint32_t end) {
	assert(isInitiated());
	for(uint32_t i = from + 1; i <= to; ++i) {
		DNAseq::value_type br = query[i-1];
		for(uint32_t j = start + 1; j <= end; ++j) {
			DNAseq::value_type bt = target[j-1];
			double s = ss->getScore(br, bt);
			double o = -ss->openGapPenalty();
			double e = -ss->extGapPenalty();
			M(i,j) = std::max({
				M(i-1,j-1) + s,
				D(i-1,j-1) + s,
				I(i-1,j-1) + s,
				0.0
			});
			I(i,j) = std::max({
				M(i-1,j) + o,
				I(i-1,j) + e
			});
			D(i,j) = std::max({
				M(i,j-1) + o,
				D(i,j-1) + e
			});
		}
	}
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_ALIGNMENTSE_H_ */
