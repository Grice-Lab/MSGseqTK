/*
 * AlignmentSE.h
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#ifndef SRC_ALIGNMENTSE_H_
#define SRC_ALIGNMENTSE_H_

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
		SeedPair(uint32_t from, uint32_t start, uint32_t len, double logP = NAN)
		: from(from), start(start), to(from + len), end(start + len), logP(logP)
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
		double logP = NAN; /* logliklihood if this SeedPair is from random match */

		/**
		 * get loglikelihood of observing this SeedPair by chance
		 */
		double loglik() const {
			return logP;
		}

		/** get the pvalue of observing this SeedPair by chance */
		double pvalue() const {
			return ::exp(loglik());
		}

		/**
		 * get in-del rate of two SeedPairs, return signed rate of the length difference devided by the longer region
		 */
		static double indelRate(const SeedPair& lhs, const SeedPair& rhs) {
			return (lhs.length() - rhs.length()) / std::max(lhs.length(), rhs.length());
		}

		/** test whether twoo SeedPairs overlap (on query or target) */
		static bool isOverlap(const SeedPair& seed1, const SeedPair& seed2) {
			return seed1.from < seed2.to && seed1.to > seed2.from ||
					seed1.start < seed2.end && seed1.end > seed2.start;
		}

		/** test whether two SeedPairs are compatitable (no overlaps) */
		static bool isCompatitable(const SeedPair& seed1, const SeedPair& seed2) {
			return !isOverlap(seed1, seed2);
		}

		/* static fields */
		static const double DEFAULT_INDEL_RATE;
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

		/**
		 * get max indel rate of between its consecutive SeedPairs
		 * @return singed max in-del rate, positive value suggesting insertion, and negative suggesting deletion
		 */
		double maxIndelRate() const;

		/** get the aggrate loglik of this SeedMatch as observing it by chance */
		double loglik() const;
	};

	typedef vector<SeedMatch> SeedMatchList;

	typedef basic_string<uint32_t> state_str;
	/* constructors */
	/** default constructor */
	AlignmentSE() = default;

	/** construct an AlignmentSE between query and target seq */
	AlignmentSE(const DNAseq& query, const DNAseq& target, const string& qname, int32_t tid, const ScoreScheme* ss)
	: query(query), target(target), qname(qname), tid(tid), ss(ss),
	  qLen(query.length()), tLen(target.length()),
	  M(qLen + 1, tLen + 1), D(qLen + 1, tLen + 1), I(qLen + 1, tLen + 1) {
		initScores();
	}

	/* member methods */
	/** test whether this alignment is initiated */
	bool isInitiated() const {
		return !query.empty() && !target.empty() && ss != nullptr;
	}

	/** initiate all score matrices */
	void initScores();

	/** calculate all scores in the entire region using Dynamic-Programming, return the alnScore as the maximum score found */
	double calculateScores() {
		calculateScores(0, qLen, 0, tLen);
		return (alnScore = M.maxCoeff(&alnTo, &alnEnd)); // determine aign 3' and score simultaneously
	}

	/** calculate all scores in the restricted "SeedMatch" regions using Dynamic-Programming, return the alnScore as the maximum score found */
	int calculateScores(const SeedMatch& seeds);

	/** backtrace alnPath */
	void backTrace();

	/** geters and setters */


	/** get align query length */
	uint32_t getAlnQLen() const {
		return alnTo - alnFrom;
	}

	/** get align target length */
	uint32_t getAlnTLen() const {
		return alnEnd - alnStart;
	}

	/** get align length, alias to getAlnTLen() */
	uint32_t getAlnLen() const {
		return getAlnTLen();
	}

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

	/**
	 * export this alignment to BAM record, given additional required fields
	 * @param tShift  shift of this alignment relative to target
	 * @param qShift  shift of this alignment relative to query
	 */
	BAM exportBAM(const string& qname, BAM::qual_str qual, int32_t tid, int32_t tShift,
			int32_t qShift = 0, uint16_t flag = 0, uint8_t mapQ = 0,
			int32_t mtid = -1, int32_t mpos = -1, int32_t isize = 0, uint64_t id = 0) const;

	/* member fields */
	DNAseq query; /* query */
	DNAseq target; /* target */
	string qname; // query name
	int32_t tid = -1;  // target id, should be determined from the database
	const ScoreScheme* ss = nullptr;

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

	/* static fileds */
	static const string STATES; // human readable STATES
	static const double DEFAULT_SEEDMATCH_EPSILON;
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


	/** get SeedMatchList from a pre-calculated MEMS and its shift of this region */
	static SeedMatchList getSeedMatchList(const MEMS& mems, uint32_t tShift,
			uint32_t qShift = 0, uint32_t maxIt = UINT16_MAX);

	/** sort SeedMatchList by their loglik, those with smaller loglik (less random) are better */
	static SeedMatchList& sortSeedMatchList(SeedMatchList& list) {
		std::sort(list.begin(), list.end(),
				[](const SeedMatch& lhs, const SeedMatch& rhs) { return lhs.loglik() < rhs.loglik(); }
		);
		return list;
	}

	/**
	 * filter SeedMatchList by using a loglik threshold that is not too much larger than the best (top) SeedMatch,
	 * assummly the input SeedMatchList is sorted
	 */
	static SeedMatchList& filterSeedMatchList(SeedMatchList& sortedList, double epsilon = DEFAULT_SEEDMATCH_EPSILON);

	/* private utility methods */
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
