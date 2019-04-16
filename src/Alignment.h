/*
 * Alignment.h
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#ifndef SRC_ALIGNMENT_H_
#define SRC_ALIGNMENT_H_

#define _USE_MATH_DEFINES
#define EIGEN_NO_DEBUG // disable eigen assertions, use one-time assertion outside loops

#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <regex> // C++11
#include "MSGseqTKConst.h"
#include "PrimarySeq.h"
#include "QualStr.h"
#include "ScoreScheme.h"
#include "SeedPair.h"
#include "SeedChain.h"
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

class Alignment; // forward declaration
class AlignmentPE; // forward declaration
typedef vector<Alignment> ALIGN_LIST;
typedef vector<AlignmentPE> PAIR_LIST;

/**
 * an alignment region of between a single-end read/query and a database/target using a seed chain
 */
class Alignment {
public:
	typedef uint32_t CIGAR_OP_TYPE;
	typedef basic_string<uint32_t> state_str;

	/* constructors */
	/** default constructor */
	Alignment() = default;

	/** construct an Alignment with all fields */
	Alignment(const PrimarySeq* read, const PrimarySeq* rcRead, const DNAseq* target,
			int32_t tid, GLoc::STRAND qStrand,
			int64_t qFrom, int64_t qTo, int64_t tStart, int64_t tEnd,
			const ScoreScheme* ss, uint8_t mapQ = INVALID_MAP_Q)
	: read(read), rcRead(rcRead), target(target),
	  tid(tid), qStrand(qStrand),
	  qFrom(qFrom), qTo(qTo), tStart(tStart), tEnd(tEnd), ss(ss), mapQ(mapQ)
	{
		assert(qFrom == 0 && qTo == read->length()); // cannot accept hard-clipped query
		init();
	}

	/** construct an Alignment between a query, database and a SeedChain */
	Alignment(const PrimarySeq* read, const PrimarySeq* rcRead, const MetaGenome& mtg,
			const ScoreScheme* ss, const SeedChain& chain, uint8_t mapQ = INVALID_MAP_Q)
	: read(read), rcRead(rcRead), tid(chain.getTid()), target(&mtg.getSeq()),
	  qStrand(chain.getStrand()),
			qFrom(0), qTo(read->length()), ss(ss), mapQ(mapQ)
	{
		init(mtg.getChromLoc(tid).getStart(), mtg.getChromLoc(tid).getEnd(), chain);
	}

	/** construct an unmapped AlignmentSE with minimum fields */
	Alignment(const PrimarySeq* read) : read(read)
	{  }

	/* member methods */
	/** getters and setters */

	int64_t getAlnEnd() const {
		return alnEnd;
	}

	int64_t getAlnFrom() const {
		return alnFrom;
	}

	const state_str& getAlnPath() const {
		return alnPath;
	}

	double getAlnScore() const {
		return alnScore;
	}

	int64_t getAlnStart() const {
		return alnStart;
	}

	int64_t getAlnTo() const {
		return alnTo;
	}

	double getLog10P() const {
		return log10P;
	}

	uint8_t getMapQ() const {
		return mapQ;
	}

	double getPostP() const {
		return postP;
	}

	int64_t getFrom() const {
		return qFrom;
	}

	GLoc::STRAND getStrand() const {
		return qStrand;
	}

	int64_t getTo() const {
		return qTo;
	}

	const PrimarySeq* getRcRead() const {
		return rcRead;
	}

	const PrimarySeq* getRead() const {
		return read;
	}

	const ScoreScheme* getSS() const {
		return ss;
	}

	const DNAseq* getTarget() const {
		return target;
	}

	int64_t getEnd() const {
		return tEnd;
	}

	int32_t getTid() const {
		return tid;
	}

	int64_t getStart() const {
		return tStart;
	}

	/** test whether this alignment is initiated */
	bool isInitiated() const {
		return read != nullptr && rcRead != nullptr && target != nullptr && ss != nullptr;
	}

	/** init Alignment with current values */
	Alignment& init();

	/** init Alignment with given MetaGenome and SeedChain */
	Alignment& init(int64_t chrStart, int64_t chrEnd, const SeedChain& chain);

	/** clear all scores to save storage (for copying/moving) */
	Alignment& clearScores();

	/**
	 * calculate all scores in a given region using Dynamic-Programming
	 * @param from, start  0-based
	 * @param to, end  1-based
	 */
	void calculateScores(int64_t from, int64_t to, int64_t start, int64_t end);

	/**
	 * calculate scores in a SeedPair region, where only M scores on the diagnal are calculated
	 */
	void calculateScores(const SeedPair& pair);

	/** calculate all scores in the entire region using Dynamic-Programming, return the alnScore as the maximum score found */
	Alignment& calculateScores() {
		calculateScores(qFrom, qTo, 0, getTLen());
		alnScore = M.maxCoeff(&alnTo, &alnEnd); // determine aign 3' and score simultaneously
		alnTo += qFrom;
		alnEnd += tStart;
		return *this;
	}

	/** calculate all scores given a SeedChain restricts using Dynamic-Programming, return the alnScore as the maximum score found */
	Alignment& calculateScores(const SeedChain& chain);

	/** backtrace alnPath */
	Alignment& backTrace();

	/** get qLen */
	int32_t getQLen() const {
		return qTo - qFrom;
	}

	/** get tLen */
	int32_t getTLen() const {
		return tEnd - tStart;
	}

	/** get 5' soft-clip size */
	int32_t getClip5Len() const {
		return alnFrom - qFrom;
	}

	/** get 3' soft-clip size */
	int32_t getClip3Len() const {
		return qTo - alnTo;
	}

	/** get all soft-clip size */
	int32_t getClipLen() const {
		return getClip5Len() + getClip3Len();
	}

	/** get align query length */
	int32_t getAlnQLen() const {
		return alnTo - alnFrom;
	}

	/** get align target length */
	int32_t getAlnTLen() const {
		return alnEnd - alnStart;
	}

	/** get align length, determined by alnPath, including M=XIDX but not H,P,N */
	int32_t getAlnLen() const;

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

	/** get MD:Z tag for this alignment */
	string getAlnMDTag() const;

	/** get overall score of this alignment, including soft-clips */
	double getScore() const {
		return alnScore - ss->getClipPenalty() * getClipLen();
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

	/** get flags based on strand information */
	uint16_t getFlag() const {
		return qStrand != GLoc::REV ? 0 : BAM_FREVERSE;
	}

	/**
	 * export core info of this alignment to a BAM record, with no aux data
	 */
	BAM exportBAM() const {
		const PrimarySeq* query = qStrand == GLoc::FWD ? read : rcRead;
		return BAM(query->getName(), getFlag(), tid, alnStart, mapQ, getAlnCigar(), getQLen(),
				dna::nt16Encode(query->getSeq()), query->getQual());
	}

private:
	/* member fields */
	const PrimarySeq* read = nullptr;
	const PrimarySeq* rcRead = nullptr;
	int32_t tid = -1;  // target id, should be determined from the database
	const DNAseq* target = nullptr;
	GLoc::STRAND qStrand = GLoc::UNK; // query strand

	int64_t qFrom = 0; // 0-based
	int64_t qTo = 0; // 1-based
	int64_t tStart = 0; // 0-based
	int64_t tEnd = 0; // 1-based

	const ScoreScheme* ss = nullptr;

	uint8_t mapQ = INVALID_MAP_Q;

	MatrixXd M; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that x[i] is aligned to y[j]
	MatrixXd D; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that x[i] is aligned to a gap (deletion for query)
	MatrixXd I; // (qLen + 1) * (tLen + 1) DP-matrix to store scores that y[j] is aligned to a gap (insertion for query)

	int64_t alnFrom = 0; // 0-based alignment from
	int64_t alnTo = 0; // 1-based alignment to
	int64_t alnStart = 0; // 0-based alignment start
	int64_t alnEnd = 0; // 1-based alignment end
	double alnScore = NAN; // alignment score
	state_str alnPath; // backtrace alignment path using cigar ops defined htslib/sam.h
//	BAM::cigar_str alnCigar; // cigar_str compressed from alnPath

	double log10P = 0;  // log10-liklihood of this alignment, given the alignment and quality
	double postP = 0;   // posterior probability of this alignment, given all candidate alignments

public:
	/* static fileds */
	static const uint8_t INVALID_MAP_Q = 0xff;
	static const uint32_t MAX_ALIGN = UINT16_MAX;
	static const double DEFAULT_INDEL_RATE;
	static const double MAX_INDEL_RATE;
	static const string STATES; // human readable STATES
	static const double DEFAULT_SCORE_REL_EPSILON; // default relative min score comparing to best score
	static const char ALIGN_GAP = '-';
	static const char ALIGN_MD_INS = '^'; // symbol use in MD:Z tag
	static const std::regex MDTAG_LEADING_PATTERN;
	static const std::regex MDTAG_MAIN_PATTERN;

	static const string MISMATCH_POSITION_TAG;
	static const string ALIGNMENT_LOG10LIK_TAG;
	static const string ALIGNMENT_POSTERIOR_PROB_TAG;

	/* static methods */
	static CIGAR_OP_TYPE matchMax(double match, double ins, double del);
	static CIGAR_OP_TYPE insMax(double match, double ins);
	static CIGAR_OP_TYPE delMax(double match, double del);

	/** get algnQLen by MD tag */
	static int32_t mdTag2alnQLen(const string& mdTag);

	/** get calculated and cleaned alignments from database and a ChainList */
	static ALIGN_LIST buildAlignments(const PrimarySeq* read, const PrimarySeq* rcRead, const MetaGenome& mtg,
			const ScoreScheme* ss, const ChainList& chains);

	/** filter candidate list of Alignments using alnScore */
	static ALIGN_LIST& filter(ALIGN_LIST& alnList, double bestFrac);

	/** evaluate each SMEM in an SMEM_LIST */
	static ALIGN_LIST& evaluate(ALIGN_LIST& alnList);

	/**
	 * calculate mapQ as posterior probability of candidate alignments using a uniform prior
	 * set the mapQ value of all alignments
	 */
	static ALIGN_LIST& calcMapQ(ALIGN_LIST& alnList);

	/** sort Alignment by mapQ/loglik decreasingly */
	static ALIGN_LIST& sort(ALIGN_LIST& alnList);

	friend class AlignmentPE;
};

/**
 * wrapper class for Alignment paired-end (PE)
 * it stores mate Alignments in pointers for performance
 * it never alter the original Alignment objects by using const pointers and keeping its own additional fields
 */
class AlignmentPE {
public:
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
		return fwdAln->tid == revAln->tid && (fwdAln->qStrand & revAln->qStrand) == 0;
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

	/** export fwd BAM */
	BAM exportFwdBAM() const {
		BAM fwdBam = fwdAln->exportBAM();
		// update paired-end info
		fwdBam.setMapQ(mapQ);
		fwdBam.setISize(isize);
		fwdBam.setFlag(getFwdFlag());
		return fwdBam;
	}

	/** export rev BAM */
	BAM exportRevBAM() const {
		BAM revBam = fwdAln->exportBAM();
		// update paired-end info
		revBam.setMapQ(mapQ);
		revBam.setISize(isize);
		revBam.setFlag(getRevFlag());
		return revBam;
	}

	/** get fwd flag */
	uint16_t getFwdFlag() const {
		return fwdAln->getFlag() | BAM_FPROPER_PAIR | BAM_FREAD1;
	}

	/** get rev flag */
	uint16_t getRevFlag() const {
		return revAln->getFlag() | BAM_FPROPER_PAIR | BAM_FREAD2;
	}

	/* member fields */
	const Alignment* fwdAln = nullptr;
	const Alignment* revAln = nullptr;
	int32_t isize = -1;

	uint8_t mapQ = Alignment::INVALID_MAP_Q;
	double postP = 0;   // posterior probability of this pair, given all candidate pairs

	/* static methods */
	/** init candidate list of Alignment pairs from fwd and rev ALIGN_LIST */
	static PAIR_LIST getPairs(const ALIGN_LIST& fwdAlnList, const ALIGN_LIST& revAlnList, uint32_t maxPair = MAX_PAIR);

	/** filter candidate pairs using pairing criteria */
	static PAIR_LIST& filter(PAIR_LIST& pairList,
			int32_t minIns, int32_t maxIns,
			bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap);

	/**
	 * calculate mapQ as posterior probability of candidate alignment pairs using a uniform prior
	 * set the mapQ value of all alignment pairs
	 */
	static PAIR_LIST& calcMapQ(PAIR_LIST& pairList);

	/** sort mapQ/loglik decreasingly for pairs */
	static PAIR_LIST& sort(PAIR_LIST& pairList);

	/** calculate insert size of two Alignment
	 * @return  distance of two alignment,
	 * or -1 if not on the same chromosome
	  */
	static int32_t getInsertSize(const Alignment& lhs, const Alignment& rhs);

	/* static fields */
	static const uint32_t MAX_PAIR = UINT16_MAX;
};

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

inline Alignment& Alignment::init() {
//	assert(isInitiated());
	M = MatrixXd::Constant(getQLen() + 1, getTLen() + 1, infV);
	I = MatrixXd::Constant(getQLen() + 1, getTLen() + 1, infV);
	D = MatrixXd::Constant(getQLen() + 1, getTLen() + 1, infV);
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

inline void Alignment::calculateScores(int64_t from, int64_t to, int64_t start, int64_t end) {
	assert(isInitiated());
	assert(qFrom <= from && to <= qTo);
	assert(tStart <= start && end <= tEnd);
	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	for(int64_t q = from; q < to; ++q) {
		int32_t i = q - qFrom + 1; // relative to score matrices
		for(int64_t t = start; t < end; ++t) {
			int32_t j = t - tStart + 1; // relative to score matrices
			double s = ss->getScore(query[q], (*target)[t]);
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
	assert(qFrom <= pair.getFrom() && pair.getTo() <= qTo && tStart <= pair.getStart() && pair.getEnd() <= tEnd);
	for(int32_t k = 0; k < pair.length(); ++k) {
		int32_t i = pair.getFrom() + k - qFrom + 1;
		int32_t j = pair.getStart() + k - tStart + 1;
		double s = ss->getMatchScore();
		if(k == 0) {
			M(i,j) = std::max({
				M(i-1,j-1) + s,
				D(i-1,j-1) + s,
				I(i-1,j-1) + s,
				0.0
			});
		}
		else
			M(i,j) = M(i-1,j-1) + s;
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
