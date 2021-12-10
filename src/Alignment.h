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
#define EIGEN_DONT_PARALLELIZE // disable eigen3 parallization, using per-read level multi-threading instead

#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <regex> // C++11
#include <utility>
#include "MSGseqTKConst.h"
#include "PrimarySeq.h"
#include "QualStr.h"
#include "ScoreScheme.h"
#include "PairingScheme.h"
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
using std::pair;

class Alignment; // forward declaration
class AlignmentPE; // forward declaration
typedef vector<Alignment> ALIGN_LIST;
typedef vector<AlignmentPE> PAIR_LIST;
typedef vector<BAM> BAM_LIST;
typedef pair<BAM, BAM> BAMPAIR;
typedef vector<BAMPAIR> BAMPAIR_LIST;

/**
 * an alignment region of between a single-end read/query and a database/target using a seed chain
 */
class Alignment {
public:
	typedef uint32_t CIGAR_OP_TYPE;
	typedef basic_string<uint32_t> state_str;

	enum MODE { GLOBAL, LOCAL };

	/* constructors */
	/** default constructor */
	Alignment() = default;

	/** construct an Alignment with all fields */
	Alignment(const PrimarySeq* read, const PrimarySeq* rcRead, const DNAseq* target,
			int32_t tid, GLoc::STRAND qStrand,
			int64_t qFrom, int64_t qTo, int64_t tStart, int64_t tEnd,
			MODE alnMode = DEFAULT_MODE, uint8_t mapQ = INVALID_MAP_Q)
	: read(read), rcRead(rcRead), tid(tid), target(target), qStrand(qStrand),
	  qFrom(qFrom), qTo(qTo), tStart(tStart), tEnd(tEnd),
	  alnMode(alnMode), mapQ(mapQ)
	{
		assert(qFrom == 0 && qTo == read->length()); // cannot accept hard-clipped query
		init();
	}

	/** construct an Alignment between a query, database and a SeedChain */
	Alignment(const PrimarySeq* read, const PrimarySeq* rcRead, const MetaGenome& mtg,
			const SeedChain& chain, MODE alnMode = DEFAULT_MODE, uint8_t mapQ = INVALID_MAP_Q)
	: read(read), rcRead(rcRead), tid(chain.getTid()), target(&mtg.getSeq(tid)), qStrand(chain.getStrand()),
	  qFrom(0), qTo(read->length()), alnMode(alnMode), mapQ(mapQ)
	{
		init(mtg.getChromLength(tid), chain);
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

	const DNAseq* getTarget() const {
		return target;
	}

	int64_t getReadLen() const {
		return read != nullptr ? read->length() : 0;
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
		return read != nullptr && rcRead != nullptr && target != nullptr;
	}

	/** init Alignment with current values */
	Alignment& init() {
		return alnMode == GLOBAL ? initNW() : initSW();
	}

	/** init Alignment for Needleman-Wunsch global alignment */
	Alignment& initNW();

	/** init Alignment for Smith-Waterman local alignment */
	Alignment& initSW();

	/** init Alignment with known chrom length and SeedChain */
	Alignment& init(int64_t chrLen, const SeedChain& chain);

	/** clear all scores to save storage (for copying/moving) */
	Alignment& clearScores();

	/**
	 * calculate all scores in a given region using Dynamic-Programming
	 * @param from, start  0-based
	 * @param to, end  1-based
	 */
	void calculateScores(int64_t from, int64_t to, int64_t start, int64_t end) {
		if(alnMode == GLOBAL)
			calculateScoresNW(from, to, start, end);
		else
			calculateScoresSW(from, to, start, end);
	}

	/**
	 * calculate all scores in a given region using Dynamic-Programming
	 * using Needleman-Wunsch global algorithm
	 * @param from, start  0-based
	 * @param to, end  1-based
	 */
	void calculateScoresNW(int64_t from, int64_t to, int64_t start, int64_t end);

	/**
	 * calculate all scores in a given region using Dynamic-Programming
	 * using Smith-Waterman local algorithm
	 * @param from, start  0-based
	 * @param to, end  1-based
	 */
	void calculateScoresSW(int64_t from, int64_t to, int64_t start, int64_t end);

	/**
	 * calculate scores in a SeedPair region, where only M scores on the diagnal are calculated
	 */
	void calculateScores(const SeedPair& pair);


	/** calculate all scores in the entire region using Dynamic-Programming, return the alnScore as the maximum score found */
	Alignment& calculateScores() {
		calculateScores(qFrom, qTo, 0, getTLen());
		return *this;
	}

	/** calculate all scores given a SeedChain restricts using Dynamic-Programming, return the alnScore as the maximum score found */
	Alignment& calculateScores(const SeedChain& chain);

	/** backtrace */
	Alignment& backTrace() {
		return alnMode == GLOBAL ? backTraceNW() : backTraceSW();
	}

	/** backtrace for Needleman-Wunsch global alignment */
	Alignment& backTraceNW();

	/** backtrace for Smith-Waterman local alignment */
	Alignment& backTraceSW();

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
		return alnScore - ss.getClipPenalty() * getClipLen();
	}

	/** evaluate this alignment */
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
	 * export core info of this alignment to a BAM record, with standard and customized aux tags set if not depend on output
	 */
	BAM exportBAM() const;

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

	MODE alnMode = DEFAULT_MODE;
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

	/* build-in alignment properties */
	double log10P = 0;  // log10-liklihood of this alignment, given the alignment and quality
	double postP = 0;   // posterior probability of this alignment, given all candidate alignments

public:
	/* static fileds */
	static ScoreScheme ss; // class-wise static score scheme object

	static const MODE DEFAULT_MODE = GLOBAL;
	static const uint8_t INVALID_MAP_Q = 0xff;
	static const uint32_t MAX_ALIGN = UINT32_MAX;
	static const double MAX_MISMATCH_RATE;
	static const double MAX_INDEL_RATE;
	static const string STATES; // human readable STATES
	static const double DEFAULT_SCORE_REL_EPSILON; // default relative min score comparing to best score
	static const char ALIGN_GAP = '-';
	static const char ALIGN_MD_INS = '^'; // symbol use in MD:Z tag
	static const std::regex MDTAG_LEADING_PATTERN;
	static const std::regex MDTAG_MAIN_PATTERN;

	/* standard aux tags */
	static const string ALIGNMENT_SCORE_TAG;
	static const string NUM_REPORTED_ALIGNMENT_TAG;
	static const string MISMATCH_POSITION_TAG;
	/* customized aux tags */
//	static const string ALIGNMENT_LENGTH_TAG;
//	static const string ALIGNMENT_INSERT_TAG;
	static const string NUM_TOTAL_ALIGNMENT_TAG;
	static const string ALIGNMENT_LOG10LIK_TAG;
	static const string ALIGNMENT_POSTERIOR_PROB_TAG;
//	static const string NUM_MISMATCHES_TAG;
//	static const string NUM_INDEL_TAG;
//	static const string ALIGNMENT_IDENTITY_TAG;

	/* static methods */
	static CIGAR_OP_TYPE matchMax(double match, double ins, double del);
	static CIGAR_OP_TYPE insMax(double match, double ins);
	static CIGAR_OP_TYPE delMax(double match, double del);

	/** get algnQLen by MD tag */
	static int32_t mdTag2alnQLen(const string& mdTag);

	/** get calculated and cleaned alignments from database and a ChainList */
	static ALIGN_LIST buildAlignments(const PrimarySeq* read, const PrimarySeq* rcRead, const MetaGenome& mtg,
			const ChainList& chains, MODE alnMode = DEFAULT_MODE);

	/** filter candidate list of Alignments using alnScore */
	static ALIGN_LIST& filter(ALIGN_LIST& alnList, double minScoreRate);

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

	/** get log10lik of this pair */
	double log10lik() const {
		return fwdAln->log10lik() + revAln->log10lik() + AlignmentPE::ps.log10lik(getInsertSize());
	}

	/** get loglik of this pair */
	double loglik() const {
		return fwdAln->loglik() + revAln->loglik() + AlignmentPE::ps.loglik(getInsertSize());
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
		fwdBam.setISize(fwdAln->alnStart < revAln->alnStart ? isize : -isize);
		fwdBam.setFlag(getFwdFlag());
		return fwdBam;
	}

	/** export rev BAM */
	BAM exportRevBAM() const {
		BAM revBam = revAln->exportBAM();
		// update paired-end info
		revBam.setMapQ(mapQ);
		revBam.setISize(fwdAln->alnStart < revAln->alnStart ? -isize : isize);
		revBam.setFlag(getRevFlag());
		return revBam;
	}

	/** get fwd flag */
	uint16_t getFwdFlag() const {
		return fwdAln->getFlag() | BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
	}

	/** get rev flag */
	uint16_t getRevFlag() const {
		return revAln->getFlag() | BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;
	}

	/* member fields */
	const Alignment* fwdAln = nullptr;
	const Alignment* revAln = nullptr;
	int32_t isize = -1;

	uint8_t mapQ = Alignment::INVALID_MAP_Q;
	double postP = 0;   // posterior probability of this pair, given all candidate pairs

	/* static methods */
	/** init candidate list of Alignment pairs from fwd and rev ALIGN_LIST */
	static PAIR_LIST getPairs(const ALIGN_LIST& fwdAlnList, const ALIGN_LIST& revAlnList);

	/** filter candidate pairs using pairing criteria */
	static PAIR_LIST& filter(PAIR_LIST& pairList,
			bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap, int64_t maxNPair = MAX_NPAIR);

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
	static PairingScheme ps; // class-wise static pairing sheme object
	static const int64_t MAX_NPAIR = 0;
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

inline Alignment& Alignment::initNW() {
	M = MatrixXd::Constant(getQLen() + 1, getTLen() + 1, infV);
	I = MatrixXd::Constant(getQLen() + 1, getTLen() + 1, infV);
	D = MatrixXd::Constant(getQLen() + 1, getTLen() + 1, infV);
	M.row(0).setZero();
	/* init first column of I */
	for(MatrixXd::Index i = 0; i < I.rows(); ++i)
		I(i, 0) = - ss.gapPenalty(i);
	return *this;
}

inline Alignment& Alignment::initSW() {
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

inline void Alignment::calculateScoresNW(int64_t from, int64_t to, int64_t start, int64_t end) {
	assert(isInitiated());
	assert(qFrom <= from && to <= qTo);
	assert(tStart <= start && end <= tEnd);
	to = std::min(to, qTo - 1);
	end = std::min(end, tEnd - 1);

	double o = -ss.openGapPenalty();
	double e = -ss.extGapPenalty();

	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	for(int64_t q = from; q <= to; ++q) {
		int32_t i = q - qFrom + 1; // relative to score matrices
		for(int64_t t = start; t <= end; ++t) {
			int32_t j = t - tStart + 1; // relative to score matrices
			double s = ss.getScore(query[q], (*target)[t]);
			M(i,j) = std::max({
				M(i-1,j-1) + s,
				D(i-1,j-1) + s,
				I(i-1,j-1) + s
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

inline void Alignment::calculateScoresSW(int64_t from, int64_t to, int64_t start, int64_t end) {
	assert(isInitiated());
	assert(qFrom <= from && to <= qTo);
	assert(tStart <= start && end <= tEnd);
	to = std::min(to, qTo - 1);
	end = std::min(end, tEnd - 1);

	double o = -ss.openGapPenalty();
	double e = -ss.extGapPenalty();

	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	for(int64_t q = from; q <= to; ++q) {
		int32_t i = q - qFrom + 1; // relative to score matrices
		for(int64_t t = start; t <= end; ++t) {
			int32_t j = t - tStart + 1; // relative to score matrices
			double s = ss.getScore(query[q], (*target)[t]);
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
//	std::cerr << "qSeq: " << read->getSeq() << std::endl;
//	std::cerr << "qFrom: " << qFrom << " qTo: " << qTo << " tStart: " << tStart << " tEnd: " << tEnd << " pair: " << pair << std::endl;
	assert(qFrom <= pair.getFrom() && pair.getTo() <= qTo && tStart <= pair.getStart() && pair.getEnd() <= tEnd);
	double o = -ss.openGapPenalty();
	double e = -ss.extGapPenalty();
	double s = ss.getMatchScore();
	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	for(int32_t k = 0; k < pair.length(); ++k) {
		int32_t i = pair.getFrom() + k - qFrom + 1;
		int32_t j = pair.getStart() + k - tStart + 1;
		if(k == 0) {
			M(i,j) = std::max({
				M(i-1,j-1) + s,
				D(i-1,j-1) + s,
				I(i-1,j-1) + s,
				0.0
			});
		}
		else {
			M(i,j) = M(i-1,j-1) + s;
		}
	}
	// process one row pass pair, if exist
	int64_t i = pair.getTo() - qFrom + 1;
	int64_t j = pair.getEnd() - tStart + 1;
	if(i <= getQLen() && j <= getTLen())
		M(i,j) = M(i-1,j-1) + s;
	if(i <= getQLen()) {
		I.row(i) = M.row(i - 1).array() + o;
	}
	// process one column pass pair, if exist
	if(j <= getTLen()) {
		D.col(j) = M.col(j - 1).array() + o;
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
