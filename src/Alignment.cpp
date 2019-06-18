/*
 * Alignment.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include "Alignment.h"

#include <cassert>
#include <climits>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>

namespace EGriceLab {
namespace MSGseqTK {
using namespace Eigen;

const double Alignment::MAX_MISMATCH_RATE = 0.2;
const double Alignment::MAX_INDEL_RATE = 0.1;
const string Alignment::STATES = BAM_CIGAR_STR;
const double Alignment::DEFAULT_SCORE_REL_EPSILON = 0.85;
const std::regex Alignment::MDTAG_LEADING_PATTERN = std::regex("^\\d+");
const std::regex Alignment::MDTAG_MAIN_PATTERN = std::regex("([A-Z]|\\^[A-Z]+)(\\d+)");
const string MISMATCH_POSITION_TAG = "MD";
const string ALIGNMENT_LOG10LIK_TAG = "XH";
const string ALIGNMENT_POSTERIOR_PROB_TAG = "XP";

Alignment& Alignment::init(int64_t chrLen, const SeedChain& chain) {
	// set up tStart and tEnd
	tStart = std::max<int64_t>(0, chain.getStart() - (1 + MAX_INDEL_RATE) * (chain.getFrom() - qFrom));
	tEnd = std::min<int64_t>(chrLen, chain.getEnd() + (1 + MAX_INDEL_RATE) * (qTo - chain.getTo()));
	init();
	return *this;
}

Alignment& Alignment::calculateScores(const SeedChain& chain) {
	assert(!chain.empty());
	/* DP at 5', if any */
	if(qFrom < chain.getFrom() && tStart < chain.getStart())
		calculateScores(qFrom, chain.getFrom(), tStart, chain.getStart());
	for(SeedChain::const_iterator seed = chain.begin(); seed < chain.end(); ++seed) {
		/* DP within this seed */
		calculateScores(*seed);
		/* DP in between this seed and next seed, if any */
		if(seed != chain.end() - 1)
			calculateScores(seed->getTo(), (seed + 1)->getFrom(), seed->getEnd(), (seed + 1)->getStart());
	}
	/* DP at 3', if any */
	if(chain.getTo() < qTo && chain.getEnd() < tEnd)
		calculateScores(chain.getTo(), qTo, chain.getEnd(), tEnd);
	return *this;
}

Alignment& Alignment::backTraceNW() {
	/* global alignment for read may exist at match or insertion state */
	alnTo = M.rows() - 1; // global alignment always exist from last column
	alnScore = infV;
	alnEnd = 0;
	uint32_t s = BAM_CMATCH; // alignment of affine DP always end and start at Match states

	/* determine the exist row and state */
	for(int64_t j = 1; j < M.cols(); ++j) {
		if(M(alnTo, j) > alnScore) {
			alnScore = M(alnTo, j);
			alnEnd = j;
			s = BAM_CMATCH;
		}
		if(I(alnTo, j) > alnScore) {
			alnScore = I(alnTo, j);
			alnEnd = j;
			s = BAM_CINS;
		}
	}

	alnTo += qFrom;
	alnEnd += tStart;
	assert(alnTo > 0 && alnEnd > 0);
	assert(alnScore != infV);
	const double o = ss->getGapOPenalty();
	const double e = ss->getGapEPenalty();

	alnPath.clear();
	alnPath.reserve(getQLen() + getTLen()); /* most time enough for back-trace */
	int64_t i = alnTo - qFrom; // 1-based
	int64_t j = alnEnd - tStart; // 1-based

	while(i > 0) {
		alnPath.push_back(s);
		switch(s) {
		case BAM_CMATCH:
			s = matchMax(M(i-1,j-1), I(i-1,j-1), D(i-1,j-1));
			i--;
			j--;
			break;
		case BAM_CINS:
			s = insMax(M(i-1,j) - o - e, I(i-1,j) - e);
			i--;
			break;
		case BAM_CDEL:
			s = delMax(M(i,j-1) - o - e, D(i,j-1) - e);
			j--;
			break;
		default:
			break;
		}
	}
	alnFrom = qFrom + i;
	alnStart = tStart + j;
	assert(!alnPath.empty() &&
			(alnPath.front() == BAM_CMATCH || alnPath.front() == BAM_CINS) &&
			(alnPath.back() == BAM_CMATCH || alnPath.back() == BAM_CINS));
	std::reverse(alnPath.begin(), alnPath.end()); // reverse alnPath
	return *this;
}

Alignment& Alignment::backTraceSW() {
	/* determine end bundaries by finding best index */
	alnScore = M.maxCoeff(&alnTo, &alnEnd); // determine aign 3' and score simultaneously
	alnTo += qFrom;
	alnEnd += tStart;
	assert(alnTo > 0 && alnEnd > 0);
	assert(alnScore != infV);

	alnPath.clear();
	alnPath.reserve(getQLen() + getTLen()); /* most time enough for back-trace */
	int64_t i = alnTo - qFrom; // 1-based
	int64_t j = alnEnd - tStart; // 1-based
	uint32_t s = BAM_CMATCH; // alignment of affine DP always end and start at Match states
	const double o = ss->getGapOPenalty();
	const double e = ss->getGapEPenalty();

	while(s != BAM_CMATCH || M(i, j) > 0) {
		alnPath.push_back(s);
		switch(s) {
		case BAM_CMATCH:
			s = matchMax(M(i-1,j-1), I(i-1,j-1), D(i-1,j-1));
			i--;
			j--;
			break;
		case BAM_CINS:
			s = insMax(M(i-1,j) - o - e, I(i-1,j) - e);
			i--;
			break;
		case BAM_CDEL:
			s = delMax(M(i,j-1) - o - e, D(i,j-1) - e);
			j--;
			break;
		default:
			break;
		}
	}
	alnFrom = qFrom + i;
	alnStart = tStart + j;
	assert(!alnPath.empty() && alnPath.front() == BAM_CMATCH && alnPath.back() == BAM_CMATCH);
	std::reverse(alnPath.begin(), alnPath.end()); // reverse alnPath
	return *this;
}

BAM::cigar_str Alignment::getAlnCigar() const {
	assert(!alnPath.empty());
	BAM::cigar_str alnCigar;
	alnCigar.reserve(alnPath.length());

	/* add 5' soft-clip, if any */
	if(getClip5Len() > 0)
		alnCigar.push_back(bam_cigar_gen(getClip5Len(), BAM_CSOFT_CLIP));

	uint32_t op = alnPath.front();
	uint32_t len = 1;
	for(uint32_t i = 1; i < alnPath.length(); ++i) {
		if(alnPath[i] != op) { /* a different op found */
			alnCigar.push_back(bam_cigar_gen(len, op)); /* add previous cigar */
			len = 1;
		}
		else
			len++;
		op = alnPath[i]; // update op
	}
	alnCigar.push_back(bam_cigar_gen(len, op)); // add last cigar op
	/* add 3' soft-clip, if any */
	if(getClip3Len() > 0)
			alnCigar.push_back(bam_cigar_gen(getClip3Len(), BAM_CSOFT_CLIP));
	return alnCigar;
}

string Alignment::getAlnMDTag() const {
	assert(!alnPath.empty());
	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	string mdTag;
	mdTag.reserve(alnPath.length());
	uint32_t prev_op = -1;
	uint32_t len = 0;
	for(int64_t k = 0, i = alnFrom, j = alnStart; k < alnPath.length(); ++k) { // i on query, j on target, k on alnPath
		nt16_t qb = query[i];
		nt16_t tb = (*target)[j];
		uint32_t op = alnPath[k];
		if(op == BAM_CMATCH) { // need distinguish EQUAL AND DIFF
			op = (qb & tb) ? BAM_CEQUAL : BAM_CDIFF; // differentiate = and X in MD tag
		}
//		fprintf(stderr, "k: %d i: %d j: %d q: %c t: %c op: %d len: %d\nquery:  %s\ntarget: %s\n", k, i, j, DNAalphabet::decode(qb), DNAalphabet::decode(tb), op, len, query->toString().c_str(), target->substr(tStart, tLen).toString().c_str());
		switch(op) {
		case BAM_CEQUAL:
			i++;
			j++;
			len++;
			break;
		case BAM_CDIFF:
			mdTag.append(boost::lexical_cast<string>(len));
			len = 0;
			mdTag.push_back(DNAalphabet::decode(tb));
			i++;
			j++;
			break;
		case BAM_CINS:
			mdTag.append(boost::lexical_cast<string>(len));
			len = 0;
			if(op != prev_op) /* open insertion */
				mdTag.push_back(ALIGN_MD_INS);
			mdTag.push_back(DNAalphabet::decode(qb));
			i++;
			break;
		case BAM_CDEL:
			j++;
			break;
		default:
			break;
		}
		prev_op = op;
	}
	/* add last length */
	mdTag.append(boost::lexical_cast<string>(len));
	return mdTag;
}

/** get decoded aligned query seq */
string Alignment::getAlnQSeq() const {
	string seq;
	seq.reserve(getAlnQLen());
	int64_t i = alnFrom; // index on query
	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	for(state_str::value_type c : alnPath) {
		switch(c) {
		case BAM_CMATCH:
			seq.push_back(dna::decode(query, i++));
			break;
		case BAM_CINS:
			seq.push_back(dna::decode(query, i++));
			break;
		case BAM_CDEL:
			seq.push_back(ALIGN_GAP);
			break;
		default:
			break;
		}
	}
	return seq;
}

/** get decoded aligned target seq */
string Alignment::getAlnTSeq() const {
	string seq;
	seq.reserve(getAlnTLen());
	int64_t j = alnStart; // index on target
	for(state_str::value_type c : alnPath) {
		switch(c) {
		case BAM_CMATCH:
			seq.push_back(dna::decode(*target, j++));
			break;
		case BAM_CINS:
			seq.push_back(ALIGN_GAP);
			break;
		case BAM_CDEL:
			seq.push_back(dna::decode(*target, j++));
			break;
		default:
			break;
		}
	}
	return seq;
}

int32_t Alignment::mdTag2alnQLen(const string& mdTag) {
	int32_t len = 0;
	std::smatch match;
	string::const_iterator searchStart = mdTag.cbegin();
	if(!std::regex_search(searchStart, mdTag.cend(), match, MDTAG_LEADING_PATTERN)) // bad formatted mdTag
		return 0;
	else
		len += boost::lexical_cast<uint32_t>(string(match[0]));
	for(searchStart = match[0].second; std::regex_search(searchStart, mdTag.cend(), match, MDTAG_MAIN_PATTERN); searchStart = match[0].second) {
		if(match.length(1) == 1) // mismatch
			len += match.length(1);
		else
			len += match.length(1) - 1;
		len += boost::lexical_cast<int32_t>(match.str(2));
	}
	return len;
}

ALIGN_LIST Alignment::buildAlignments(const PrimarySeq* read, const PrimarySeq* rcRead, const MetaGenome& mtg,
		const ScoreScheme* ss, const ChainList& chains, MODE alnMode) {
	ALIGN_LIST alnList;
	alnList.reserve(chains.size());
	for(const SeedChain& chain : chains)
		alnList.push_back(Alignment(read, rcRead, mtg, ss, chain, alnMode).calculateScores(chain).backTrace().clearScores());
	return alnList;
}

int32_t Alignment::getAlnLen() const {
	uint32_t alnLen = 0;
	for(state_str::value_type s : alnPath) {
		if(bam_cigar_type(s))
			alnLen++;
	}
	return alnLen;
}

Alignment& Alignment::evaluate() {
	if(alnPath.empty())
		*this;
	log10P = 0;
	const DNAseq& query = qStrand == GLoc::FWD ? read->getSeq() : rcRead->getSeq();
	const QualStr& qual = qStrand == GLoc::FWD ? read->getQual() : rcRead->getQual();
	/* process 5' soft-clips, if any */
	for(int64_t i = qFrom; i < alnFrom; ++i)
		log10P += - ss->getClipPenalty();

	for(int64_t k = 0, i = alnFrom, j = alnStart; k < alnPath.length(); ++k) { /* k on alnPath, i on query, j on target */
		state_str::value_type s = alnPath[k];
//		fprintf(stderr, "k: %d i: %d j: %d s: %c q: %c t: %c qual: %d log10P: %g\n", k, i, j, bam_cigar_opchr(s), DNAalphabet::decode(query[i]), DNAalphabet::decode((*target)[j]), qual[i], log10P);
		switch(s) {
		case BAM_CMATCH:
			if(query[i] & (*target)[j]) // IUPAC match (BAM_CEQUAL)
				log10P += ::log10(1 - quality::phredQ2P(qual[i]));
			else
				log10P += qual[i] / quality::PHRED_SCALE;
			i++;
			j++;
			break;
		case BAM_CINS:
			log10P += qual[i] / quality::PHRED_SCALE; // basic mismatch penalty
			if(k == 0 || alnPath[k - 1] != BAM_CINS) // gap open
				log10P += - ss->getGapOPenalty(); // additional penalty
			log10P += -ss->getGapEPenalty();
			i++;
			break;
		case BAM_CDEL:
			log10P += qual[i] / quality::PHRED_SCALE; // basic mismatch penalty
			if(k == 0 || alnPath[k - 1] != BAM_CDEL) // gap open
				log10P += -ss->getGapOPenalty();
			log10P += -ss->getGapEPenalty();
			j++;
			break;
		default:
			break;
		}
	}

	/* process 3' soft-clips, if any */
	for(int64_t i = alnTo; i < qTo; ++i)
		log10P += - ss->getClipPenalty();
	return *this;
}

ALIGN_LIST& Alignment::calcMapQ(ALIGN_LIST& alnList) {
	if(alnList.empty())
		return alnList;
	const size_t N = alnList.size();
	VectorXd pr(N);    // alignment probability
	for(size_t i = 0; i < N; ++i)
		pr(i) = alnList[i].loglik();
	double maxV = pr.maxCoeff();
	double scale = maxV != infV && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;
	pr = (pr.array() + scale).exp();
	assert(pr.sum() > 0);
	VectorXd postP = pr / pr.sum();
	/* assign postP and mapQ */
	for(size_t i = 0; i < N; ++i) {
		alnList[i].postP = postP(i);
		double mapQ = quality::phredP2Q(1 - postP(i));
		assert(!std::isnan(mapQ));
		alnList[i].mapQ = std::min(::round(mapQ), static_cast<double> (quality::MAX_Q_SCORE));
	}
	return alnList;
}

PAIR_LIST& AlignmentPE::calcMapQ(PAIR_LIST& pairList) {
	if(pairList.empty())
		return pairList;
	const size_t N = pairList.size();
//	VectorXd prior = VectorXd::Ones(N); // use a uniform prior
	VectorXd pr(N);    // alignment probability
	for(size_t i = 0; i < N; ++i)
		pr(i) = pairList[i].loglik();
	double maxV = pr.maxCoeff();
	double scale = maxV != infV && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;
	pr = (pr.array() + scale).exp();
	/* get postP */
	VectorXd postP(pr.rows());
	if(pr.sum() > 0)
		postP = pr / pr.sum(); // uniform prior ignored
	else
		postP.fill(1.0 / pr.rows()); // all equal zeros after scale
	/* assign postP and mapQ */
	for(size_t i = 0; i < N; ++i) {
		pairList[i].postP = postP(i);
		double mapQ = quality::phredP2Q(1 - postP(i));
		assert(!std::isnan(mapQ));
		pairList[i].mapQ = std::min(::round(mapQ), static_cast<double> (quality::MAX_Q_SCORE));
	}
	return pairList;
}

int32_t AlignmentPE::getInsertSize(const Alignment& lhs, const Alignment& rhs) {
	if(lhs.tid != rhs.tid)
		return -1;
	else
		return std::max(lhs.alnEnd, rhs.alnEnd) - std::min(lhs.alnStart, rhs.alnStart);
}

ALIGN_LIST& Alignment::filter(ALIGN_LIST& alnList, double bestFrac) {
	if(alnList.size() <= 1)
		return alnList;
	/* find best alignment by alnScore */
	ALIGN_LIST::const_iterator bestAln = std::max_element(alnList.cbegin(), alnList.cend(), [] (const Alignment& lhs, const Alignment& rhs) { return lhs.alnScore > rhs.alnScore; });
	/* filter alignment by alnScore */
	alnList.erase(std::remove_if(alnList.begin(), alnList.end(), [&] (const Alignment& aln) { return aln.alnScore < bestAln->alnScore * bestFrac; }), alnList.end());
	return alnList;
}

PAIR_LIST& AlignmentPE::filter(PAIR_LIST& pairList,
		int32_t minIns, int32_t maxIns,
		bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap) {
	if(!(minIns == 0 && maxIns == 0)) // filter by insert size
		pairList.erase(std::remove_if(pairList.begin(), pairList.end(), [=] (const AlignmentPE& pair) { return !(minIns <= pair.getInsertSize() && pair.getInsertSize() <= maxIns); }), pairList.end());
	if(noDiscordant) // filter discordant pairs
		pairList.erase(std::remove_if(pairList.begin(), pairList.end(), [] (const AlignmentPE& pair) { return !pair.isConcordant(); }), pairList.end());
	if(noTailOver) // filter tail-overlap pairs
		pairList.erase(std::remove_if(pairList.begin(), pairList.end(), [] (const AlignmentPE& pair) { return pair.isTailOver(); }), pairList.end());
	if(noContain) // filter containing pairs
		pairList.erase(std::remove_if(pairList.begin(), pairList.end(), [] (const AlignmentPE& pair) { return pair.isContained(); }), pairList.end());
	if(noOverlap) // filter overlap pairs
		pairList.erase(std::remove_if(pairList.begin(), pairList.end(), [] (const AlignmentPE& pair) { return pair.isOverlap(); }), pairList.end());
	return pairList;
}

PAIR_LIST AlignmentPE::getPairs(const ALIGN_LIST& fwdAlnList, const ALIGN_LIST& revAlnList, uint32_t maxPair) {
	PAIR_LIST pairList;
	pairList.reserve(std::min<size_t>(fwdAlnList.size() * revAlnList.size(), maxPair));
	for(ALIGN_LIST::const_iterator fwdAln = fwdAlnList.begin(); fwdAln < fwdAlnList.end() && pairList.size() < maxPair; ++fwdAln) {
		for(ALIGN_LIST::const_iterator revAln = revAlnList.begin(); revAln < revAlnList.end() && pairList.size() < maxPair; ++revAln) {
			if((fwdAln->qStrand & revAln->qStrand) == 0) // AlignmentPE must be on different strand
				pairList.push_back(AlignmentPE(&*fwdAln, &*revAln));
		}
	}
	return pairList;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

