/*
 * Alignment.h
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#ifndef SRC_ALIGNMENT_H_
#define SRC_ALIGNMENT_H_

#define _USE_MATH_DEFINES

#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
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

struct Alignment; // forward declaration
struct AlignmentPE; // forward declaration
typedef vector<Alignment> ALIGN_LIST;
typedef vector<AlignmentPE> PAIR_LIST;

/**
 * an alignment region of between a single-end read/query and a database/target
 */
struct Alignment {
	typedef uint32_t CIGAR_OP_TYPE;
	/* nested types and enums */
	/** A SeedPair is a pair of exact matched seed region between query and target */
	struct SeedPair {
		/* constructors */
		/** default constructor */
		SeedPair() = default;

		/** construct from given values */
		SeedPair(uint64_t from, uint64_t start, uint32_t len, int32_t tid, double logP)
		: from(from), to(from + len), start(start), end(start + len), tid(tid), logP(logP)
		{ 	}

		/* member methods */
		uint32_t length() const {
			return end - start;
		}

		/** test whether this SeedPair is valid */
		bool isValid() const {
			return tid >= 0 && !std::isnan(logP);
		}

		/* member fields */
		uint64_t from = 0; /* 0-based query position */
		uint64_t to = 0;   /* 1-based query position */
		uint64_t start = 0; /* 0-based target position */
		uint64_t end = 0; /* 1-based target position */
		int32_t tid = -1; /* tid as chromId/locId determined by database */
		double logP = NAN; /* log-probability (loglik) of observing this SeedPair by chance */

		/** test whether two SeedPairs overlap (on target) */
		static bool isOverlap(const SeedPair& lhs, const SeedPair& rhs) {
			return lhs.tid == rhs.tid && lhs.start < rhs.end && lhs.end > rhs.start;
		}

		/** get total indels between two seeds, positive values indicates insertions, and negative values gives deletions */
		static int64_t numGaps(const SeedPair& lhs, const SeedPair& rhs) {
			return (rhs.from - lhs.to) - (rhs.start - lhs.end);
		}

		/** test whether two SeedPairs are compatitable */
		static bool isCompatitable(const SeedPair& lhs, const SeedPair& rhs) {
			return lhs.tid == rhs.tid && lhs.to < rhs.from && lhs.end < rhs.start;
		}

		/** test whether two SeedPairs are compatitable and ordered and with not too much indels */
		static bool isCompatitable(const SeedPair& lhs, const SeedPair& rhs, uint64_t maxIndel) {
			return isCompatitable(lhs, rhs) &&
					::abs(numGaps(lhs, rhs)) <= maxIndel;
		}
	};

	/**
	 * a SeedMatch is an ordered lists of SeedPairs
	 */
	struct SeedMatch : public std::vector<SeedPair> {
		/* member methods */
		/** test whether this SeedMatch is compatitable */
		bool isCompatitable() const;

		/** test whether this SeedMatch is compatitable */
		bool isCompatitable(uint64_t maxIndel) const;

		/** get total length of SeedPairs within */
		uint32_t length() const;

		/** get the from position as the first from */
		uint64_t getFrom() const {
			return front().from;
		}

		/** get the to position as the last to */
		uint64_t getTo() const {
			return back().to;
		}

		/** get start position of this SeedMatch as the first SeedPair start */
		uint64_t getStart() const {
			return front().start;
		}

		/** get end position of this SeedMatch as the last SeedPair end */
		uint64_t getEnd() const {
			return back().end;
		}

		/** get tid of this SeedMatch */
		int32_t getTId() const {
			return front().tid;
		}

		/** get loglik of this SeedMatch */
		double loglik() const {
			double logP = 0;
			for(const SeedPair& seed : *this)
				logP += seed.logP;
			return logP;
		}

		/** filter this SeedMatch by removing invalid seeds */
		SeedMatch& filter();

		/** filter this SeedMatch using greedy algorithm,
		 * by progressively remove bad SeedPair from it */
		SeedMatch& filter(uint64_t maxIndel);
	};

	typedef vector<SeedMatch> SeedMatchList;
	typedef basic_string<uint32_t> state_str;

	/* constructors */
	/** default constructor */
	Alignment() = default;

	/** construct an AlignmentSE with all fields */
	Alignment(const DNAseq* query, const DNAseq* target, const QualStr* qual, const string* qname, int32_t tid,
			uint64_t qFrom, uint64_t qTo, uint64_t tStart, uint64_t tEnd,
			const ScoreScheme* ss, uint16_t flag, uint8_t mapQ, int32_t mtid, int32_t mpos, int32_t isize, uint32_t id)
	: query(query), target(target), qual(qual), qname(qname), tid(tid),
	  qFrom(qFrom), qTo(qTo), tStart(tStart), tEnd(tEnd), qLen(qTo - qFrom), tLen(tEnd - tStart),
	  ss(ss), flag(flag), mapQ(mapQ), mtid(mtid), mpos(mpos), isize(isize), id(id),
	  M(MatrixXd::Constant(qLen + 1, tLen + 1, infV)),
	  I(MatrixXd::Constant(qLen + 1, tLen + 1, infV)),
	  D(MatrixXd::Constant(qLen + 1, tLen + 1, infV)) {
		assert(qFrom == 0 && qTo == query->length()); // cannot accept hard-clipped query
		initScores();
	}

	/** construct an AlignmentSE with required fields */
	Alignment(const DNAseq* query, const DNAseq* target, const QualStr* qual, const string* qname, int32_t tid,
			uint64_t qFrom, uint64_t qTo, uint64_t tStart, uint64_t tEnd,
			const ScoreScheme* ss, uint16_t flag)
	: query(query), target(target), qual(qual), qname(qname), tid(tid),
	  qFrom(qFrom), qTo(qTo), tStart(tStart), tEnd(tEnd), qLen(qTo - qFrom), tLen(tEnd - tStart),
	  ss(ss), flag(flag),
	  M(MatrixXd::Constant(qLen + 1, tLen + 1, infV)),
	  I(MatrixXd::Constant(qLen + 1, tLen + 1, infV)),
	  D(MatrixXd::Constant(qLen + 1, tLen + 1, infV))
	{
		assert(qFrom == 0 && qTo == query->length()); // cannot accept hard-clipped query
		initScores();
	}

	/* member methods */
	/** test whether this alignment is initiated */
	bool isInitiated() const {
		return query != nullptr && target != nullptr && qual != nullptr && ss != nullptr;
	}

	/** initiate all score matrices */
	Alignment& initScores();

	/** clear all scores to save storage (for copying/moving) */
	Alignment& clearScores();

	/**
	 * calculate all scores in a given region using Dynamic-Programming
	 * @param from, start  0-based
	 * @param to, end  1-based
	 */
	void calculateScores(uint64_t from, uint64_t to, uint64_t start, uint64_t end);

	/**
	 * calculate scores in a SeedPair region, where only M scores on the diagnal are calculated
	 */
	void calculateScores(const SeedPair& pair);

	/** calculate all scores in the entire region using Dynamic-Programming, return the alnScore as the maximum score found */
	Alignment& calculateScores() {
		calculateScores(0, qLen, 0, tLen);
		alnScore = M.maxCoeff(&alnTo, &alnEnd); // determine aign 3' and score simultaneously
		alnTo += qFrom;
		alnEnd += tStart;
		return *this;
	}

	/** calculate all scores in the restricted "SeedMatch" regions using Dynamic-Programming, return the alnScore as the maximum score found */
	Alignment& calculateScores(const SeedMatch& seeds);

	/** backtrace alnPath */
	Alignment& backTrace();

	/** get strand of this Alignment */
	MEM::STRAND getStrand() const {
		return flag & BAM_FREVERSE == 0 ? MEM::FWD : MEM::REV;
	}

	/** get 5' soft-clip size */
	uint32_t getClip5Len() const {
		return alnFrom - qFrom;
	}

	/** get 3' soft-clip size */
	uint32_t getClip3Len() const {
		return qTo - alnTo;
	}

	/** get all soft-clip size */
	uint32_t getClipLen() const {
		return getClip5Len() + getClip3Len();
	}

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

	/** get overall score of this alignment, including soft-clips */
	double getScore() const {
		return alnScore - ss->clipPenalty * getClipLen();
	}

	/** evaluate this alignment log-liklihood using seq, align-path and quality */
	Alignment& evaluate();

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

	/** set flag */
	Alignment& setFlag(uint16_t flag) {
		this->flag |= flag;
		return *this;
	}

	/**
	 * export core info of this alignment to a BAM record, with no aux data
	 */
	BAM exportBAM() const {
		return BAM(*qname, flag, tid, alnStart, mapQ,
				getAlnCigar(), qLen, nt16Encode(*query), *qual,
				mtid, mpos, isize, id);
	}

	/* member fields */
	/* htslib required fields */
	const DNAseq* query = nullptr;
	const DNAseq* target = nullptr;
	const QualStr* qual = nullptr; // query qual
	const string* qname = nullptr; // query name
	int32_t tid = -1;  // target id, should be determined from the database

	uint64_t qFrom = 0; // 0-based
	uint64_t qTo = 0; // 1-based
	uint64_t tStart = 0; // 0-based
	uint64_t tEnd = 0; // 1-based
	uint32_t qLen = 0; // query length (htslib query length is restricted to max chrom size (UITN32_MAX)
	uint32_t tLen = 0; // target length

	const ScoreScheme* ss = nullptr;

	uint16_t flag = 0;
	uint8_t mapQ = INVALID_MAP_Q;
	int32_t mtid = -1;
	int32_t mpos = -1;
	int32_t isize = 0;
	uint64_t id = 0;

	MatrixXd M; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that x[i] is aligned to y[j]
	MatrixXd D; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that x[i] is aligned to a gap (deletion for query)
	MatrixXd I; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that y[j] is aligned to a gap (insertion for query)

	uint64_t alnFrom = 0; // 0-based alignment from
	uint64_t alnTo = 0; // 1-based alignment to
	uint64_t alnStart = 0; // 0-based alignment start
	uint64_t alnEnd = 0; // 1-based alignment end
	double alnScore = NAN; // alignment score
	state_str alnPath; // backtrace alignment path using cigar ops defined htslib/sam.h
//	BAM::cigar_str alnCigar; // cigar_str compressed from alnPath

	double log10P = 0;  // log10-liklihood of this alignment, given the alignment and quality
	double postP = 0;   // posterior probability of this alignment, given all candidate alignments

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

	/** get candidate Alignments from SeedMatch list */
	static ALIGN_LIST getAlignments(const ScoreScheme& ss, const MetaGenome& mtg, const PrimarySeq& read,
			const SeedMatchList& seedMatches, MEM::STRAND strand);

	/** filter candidate list of Alignments using alnScore */
	static ALIGN_LIST& filter(ALIGN_LIST& alnList, double bestFrac);

	static ALIGN_LIST& evaluate(ALIGN_LIST& alnList);

	/**
	 * calculate mapQ as posterior probability of candidate alignments using a uniform prior
	 * set the mapQ value of all alignments
	 */
	static ALIGN_LIST& calcMapQ(ALIGN_LIST& alnList);

	/** sort Alignment by mapQ/loglik decreasingly */
	static ALIGN_LIST& sort(ALIGN_LIST& alnList);
};

/**
 * wrapper class for Alignment paired-end (PE)
 * it stores mate Alignments in pointers for performance
 * it never alter the original Alignment objects by using const pointers and keeping its own additional fields
 */
struct AlignmentPE {
	/* constructors */
	/** default constructor */
	AlignmentPE() = default;

	/** construct AlignmentPE given both mates */
	AlignmentPE(const Alignment* fwdAln, const Alignment* revAln) : fwdAln(fwdAln), revAln(revAln)
	{
		setInsertSize(getInsertSize(*fwdAln, *revAln));
	}

	/* member methods */
	/** test whether this AlignmentPE is concordance */
	bool isConcordant() const {
		return fwdAln->tid == revAln->tid && fwdAln->getStrand() != revAln->getStrand();
	}

	/** test whethther two mates is tail-overlap through each other */
	bool isTailOver() const {
		return fwdAln->alnStart <= revAln->alnStart && fwdAln->alnEnd > revAln->alnStart || /* fwd-rev */
			   revAln->alnStart <= fwdAln->alnStart && revAln->alnEnd > fwdAln->alnStart; /* rev-fwd */
	}

	/** test whether one mate is contained in the other */
	bool isContained() const {
		return fwdAln->alnStart <= revAln->alnStart && fwdAln->alnEnd >= revAln->alnEnd || /* fwd contains rev */
			   revAln->alnStart <= fwdAln->alnStart && revAln->alnEnd >= fwdAln->alnEnd;
	}

	/** test whether two mates are overlapping to each other */
	bool isOverlap() const {
		return fwdAln->alnStart < revAln->alnEnd && fwdAln->alnEnd > revAln->alnStart;
	}

	/** get alnScore of this pair */
	double getAlnScore() const {
		return fwdAln->alnScore + revAln->alnScore;
	}

	/** get overal score of this pair */
	double getScore() const {
		return fwdAln->getScore() + revAln->getScore();
	}

	double log10lik() const {
		return fwdAln->log10lik() + revAln->log10lik();
	}

	/** get loglik of this pair */
	double loglik() const {
		return fwdAln->loglik() + revAln->loglik();
	}

	/** get insert size */
	int32_t getInsertSize() const {
		return isize;
	}

	/** set insert size */
	void setInsertSize(int32_t isize) {
		this->isize = isize;
	}

	/* member fields */
	const Alignment* fwdAln = nullptr;
	const Alignment* revAln = nullptr;
	int32_t isize = -1;

	uint8_t mapQ = Alignment::INVALID_MAP_Q;
	double postP = 0;   // posterior probability of this pair, given all candidate pairs

	/* static methods */
	/** init candidate list of Alignment pairs from fwd and rev ALIGN_LIST */
	static PAIR_LIST getPairs(ALIGN_LIST& fwdAlnList, ALIGN_LIST& revAlnList);

	/** filter candidate pairs using pairing criteria */
	static PAIR_LIST& filter(PAIR_LIST& pairList,
			int32_t minIns, int32_t maxIns,
			bool noMixed, bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap);

	/**
	 * calculate mapQ as posterior probability of candidate alignment pairs using a uniform prior
	 * set the mapQ value of all alignment pairs
	 */
	static PAIR_LIST& calcMapQ(PAIR_LIST& pairList);

	/** sort mapQ/loglik decreasingly for pairs */
	static PAIR_LIST& sort(PAIR_LIST& pairList);

	/** calculate insert size of two Alignment
	 * @return  zero if overlap,
	 * or -1 if not on the same chromosome/target
	  */
	static int32_t getInsertSize(const Alignment& lhs, const Alignment& rhs);
};

inline Alignment::SeedMatch& Alignment::SeedMatch::filter() {
	erase(std::remove_if(begin(), end(), [] (const SeedPair& seed) { return !seed.isValid(); }), end());
	return *this;
}

inline Alignment::CIGAR_OP_TYPE Alignment::matchMax(double match, double ins, double del) {
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

inline Alignment::CIGAR_OP_TYPE Alignment::insMax(double match, double ins) {
	return match > ins ? BAM_CMATCH : BAM_CINS;
}

inline Alignment::CIGAR_OP_TYPE Alignment::delMax(double match, double del) {
	return match > del ? BAM_CMATCH : BAM_CDEL;
}

inline Alignment& Alignment::initScores() {
//	assert(isInitiated());
	M.row(0).setZero();
	M.col(0).setZero();
	return *this;
}

inline Alignment& Alignment::clearScores() {
	M = MatrixXd();
	I = MatrixXd();
	D = MatrixXd();
	return *this;
}

inline void Alignment::calculateScores(uint64_t from, uint64_t to, uint64_t start, uint64_t end) {
	assert(isInitiated());
	for(uint64_t q = from; q < to; ++q) {
		uint32_t i = q - qFrom + 1; // relative to score matrices
		for(uint64_t t = start; t < end; ++t) {
			uint32_t j = t - tStart + 1; // relative to score matrices
			double s = ss->getScore((*query)[q], (*target)[t]);
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

inline void Alignment::calculateScores(const SeedPair& pair) {
	assert(isInitiated());
	for(uint32_t k = 0; k < pair.length(); ++k) {
		uint64_t q = pair.from + k;
		uint64_t t = pair.start + k;
		uint32_t i = q - qFrom + 1;
		uint32_t j = t - tStart + 1;
		M(i,j) = M(i-1,j-1) + ss->getScore((*query)[q], (*target)[t]);
	}
}

inline ALIGN_LIST& Alignment::evaluate(ALIGN_LIST& alnList) {
	/* evaluate filtered alignments */
	for(Alignment& aln : alnList)
		aln.evaluate();
	return alnList;
}

inline ALIGN_LIST& Alignment::sort(ALIGN_LIST& alnList) {
	/* sort candidates by their log10-liklihood descreasingly (same order as postP with uniform prior */
	std::sort(alnList.begin(), alnList.end(), [] (const Alignment& lhs, const Alignment& rhs) { return lhs.log10lik() > rhs.log10lik(); });
	return alnList;
}

inline PAIR_LIST& AlignmentPE::sort(PAIR_LIST& pairList) {
	/* sort candidates by their log10-liklihood descreasingly (same order as postP with uniform prior */
	std::sort(pairList.begin(), pairList.end(), [] (const AlignmentPE& lhs, const AlignmentPE& rhs) { return lhs.log10lik() > rhs.log10lik(); });
	return pairList;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_ALIGNMENT_H_ */
