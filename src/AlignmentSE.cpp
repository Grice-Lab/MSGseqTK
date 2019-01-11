/*
 * AlignmentSE.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include <cassert>
#include <climits>
#include <algorithm>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include "AlignmentSE.h"

namespace EGriceLab {
namespace MSGseqTK {
using namespace Eigen;

const double AlignmentSE::DEFAULT_INDEL_RATE = 0.02;
const double AlignmentSE::MAX_INDEL_RATE = 0.15;
const string AlignmentSE::STATES = BAM_CIGAR_STR;
const double AlignmentSE::DEFAULT_SCORE_REL_EPSILON = 0.85;
const boost::regex AlignmentSE::MDTAG_LEADING_PATTERN = boost::regex("^\\d+");
const boost::regex AlignmentSE::MDTAG_MAIN_PATTERN = boost::regex("([A-Z]|\\^[A-Z]+)(\\d+)");

AlignmentSE& AlignmentSE::calculateScores(const SeedMatch& seeds) {
	assert(!seeds.empty());
//	assert(seeds.isCompatitable());
	/* DP at 5', if any */
	calculateScores(qFrom, seeds.front().from, tStart, seeds.front().start);

	for(SeedMatch::const_iterator seed = seeds.begin(); seed < seeds.end(); ++seed) {
		/* DP within this seed */
		calculateScores(*seed);
		/* DP in between this seed and next seed, if any */
		if(seed != seeds.end() - 1)
			calculateScores(seed->to, (seed + 1)->from, seed->end, (seed + 1)->start);
	}
	/* DP at 3', if any */
	calculateScores(seeds.back().to, qTo, seeds.back().end, tEnd); /* could be empty loop */
	alnScore = M.maxCoeff(&alnTo, &alnEnd); // determine aign 3' and score simultaneously
	alnTo += qFrom;
	alnEnd += tStart;

	return *this;
}

AlignmentSE& AlignmentSE::backTrace() {
	assert(alnScore != infV);
//	cerr << "qname: " << qname << " tid: " << tid << " alnFrom: " << alnFrom << " alnTo: " << alnTo << " alnStart: " << alnStart << " alnEnd: " << alnEnd << endl;
	assert(alnTo > 0 && alnEnd > 0);
	alnPath.clear();
	alnPath.reserve(qLen + tLen); /* most time enough for back-trace */
	int64_t i = alnTo - qFrom; // 1-based
	int64_t j = alnEnd - tStart; // 1-based
	uint32_t s = BAM_CMATCH; // local alignment of affine DP always end and start at Match states

	while(i >= 0 && j >= 0) {
		if(s == BAM_CMATCH && M(i,j) == 0)
			break;
		alnPath.push_back(s);
		switch(s) {
		case BAM_CMATCH:
			s = matchMax(M(i-1,j-1), I(i-1,j-1), D(i-1,j-1));
			i--;
			j--;
			break;
		case BAM_CINS:
			s = insMax(M(i-1,j) - ss->openGapPenalty(), I(i-1,j) - ss->extGapPenalty());
			i--;
			break;
		case BAM_CDEL:
			s = delMax(M(i,j-1) - ss->openGapPenalty(), D(i,j-1) - ss->extGapPenalty());
			j--;
			break;
		default:
			break;
		}
	}
	alnFrom = qFrom + i;
	alnStart = tStart + j;
	assert(alnPath.front() == BAM_CMATCH && alnPath.back() == BAM_CMATCH);
	std::reverse(alnPath.begin(), alnPath.end()); // reverse alnPath
//	if(!(alnPath.front() == BAM_CMATCH && alnPath.back() == BAM_CMATCH)) {
//		fprintf(stderr, "alnFrom: %d alnTo: %d alnStart: %d alnEnd: %d alnScore %g alnCigar: %s\nquery:  %s\ntarget: %s\n",
//				alnFrom, alnTo, alnStart, alnEnd, alnScore, BAM::decodeCigar(getAlnCigar()).c_str(),
//				query->toString().c_str(), target->toString().c_str());
//		std::cerr << "M: " << std::endl << M << std::endl;
//	}
	return *this;
}

BAM::cigar_str AlignmentSE::getAlnCigar() const {
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

string AlignmentSE::getAlnMDTag() const {
	assert(!alnPath.empty());
	string mdTag;
	mdTag.reserve(alnPath.length());
	uint32_t prev_op = -1;
	uint32_t len = 0;
	for(uint64_t k = 0, i = alnFrom, j = alnStart; k < alnPath.length(); ++k) { // i on query, j on target, k on alnPath
		nt16_t qb = (*query)[i];
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

BAM::seq_str AlignmentSE::nt16Encode(const DNAseq& seq) {
	const uint32_t L = seq.length();
	BAM::seq_str seq16((L + 1) / 2, 0); // ceil(L / 2)
	for(uint32_t i = 0; i < L; i += 2)
		seq16[i / 2] = BAM::nt16Encode(seq.decode(i)) << 4 | BAM::nt16Encode(seq.decode(i+1)); // seq[i+1] always valid with the null terminal
	return seq16;
}

bool AlignmentSE::SeedMatch::isCompatitable() const {
	if(size() <= 1) // less than 1 SeedPair
		return true;

	for(SeedMatch::const_iterator seed = begin(); seed < end() - 1; ++seed)
		if(!SeedPair::isCompatitable(*seed, *(seed + 1)))
			return false;
	return true;
}

bool AlignmentSE::SeedMatch::isCompatitable(uint64_t maxIndel) const {
	std::cerr << "testing seedmatch of size " << size() << std::endl;
	for(const SeedPair& seed : *this) {
		fprintf(stderr, "from: %d to: %d start: %d end: %d tid: %d\n", seed.from, seed.to, seed.start, seed.end, seed.tid);
	}
	if(size() <= 1) // less than 1 SeedPair
		return true;

	for(SeedMatch::const_iterator seed = begin(); seed < end() - 1; ++seed)
		if(!SeedPair::isCompatitable(*seed, *(seed + 1), maxIndel))
			return false;
	return true;
}

uint32_t AlignmentSE::SeedMatch::length() const {
	uint32_t len = 0;
	for(const SeedMatch::value_type& seed : *this)
		len += seed.length();
	return len;
}

AlignmentSE::SeedMatchList AlignmentSE::getSeedMatchList(const MetaGenome& mtg, const MEMS& mems,
		uint32_t maxIt) {
	/* get a raw SeedMatchList to store all matches of each MEMS */
	const size_t N = mems.size();
	SeedMatchList rawList, outputList;
	rawList.resize(N);
	for(size_t i = 0; i < N; ++i) {
		const MEM& mem = mems[i];
		for(const Loc& loc : mem.locs) {
			int32_t tid = mtg.getChromIndex(loc.start); // use chrom index as tid
			assert(loc.start - mtg.getChromStart(tid) <= UINT32_MAX); // BAM file only support up-to UINT32_MAX chrom size
			/* rawList SeedPair start/end is relative to the chromosome */
			rawList[i].push_back(SeedPair(mem.from, loc.start, mem.length(), tid));
		}
	}

	vector<size_t> idx(N, 0); // index to keep track of next element in each of the N SeedMatch
	/* non-recursive algorithm to get SeedMatchList by randomly picking up elements from rawList */
	while(outputList.size() <= maxIt) {
		/* get current combination */
		SeedMatch combination;
		combination.reserve(N);
		for(size_t i = 0; i < N; ++i)
			combination.push_back(rawList[i][idx[i]]);
		if(combination.isCompatitable(mems.getSeq()->length() * MAX_INDEL_RATE))
			outputList.push_back(combination); // add this combination

		/* find the rightmost SeedMatch that has more elemtns left after current element */
		int64_t next = N - 1;
		while(next >= 0 && idx[next] + 1 >= rawList[next].size())
			next--;
		if(next < 0) // so such SeedMatch found, all combination explored
			break;
		idx[next]++;  // if found move to next element in that array

		/* reset index right of next */
		for(size_t i = next + 1; i < N; ++i)
			idx[i] = 0;
	}

	return outputList;
}

/** get decoded aligned query seq */
string AlignmentSE::getAlnQSeq() const {
	string seq;
	seq.reserve(getAlnQLen());
	uint64_t i = alnFrom; // index on query
	for(state_str::value_type c : alnPath) {
		switch(c) {
		case BAM_CMATCH:
			seq.push_back(query->decode(i++));
			break;
		case BAM_CINS:
			seq.push_back(query->decode(i++));
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
string AlignmentSE::getAlnTSeq() const {
	string seq;
	seq.reserve(getAlnTLen());
	uint64_t j = alnStart; // index on target
	for(state_str::value_type c : alnPath) {
		switch(c) {
		case BAM_CMATCH:
			seq.push_back(target->decode(j++));
			break;
		case BAM_CINS:
			seq.push_back(ALIGN_GAP);
			break;
		case BAM_CDEL:
			seq.push_back(target->decode(j++));
			break;
		default:
			break;
		}
	}
	return seq;
}

uint32_t AlignmentSE::mdTag2alnQLen(const string& mdTag) {
	uint32_t len = 0;
	boost::smatch match;
	string::const_iterator searchStart = mdTag.cbegin();
	if(!boost::regex_search(searchStart, mdTag.cend(), match, MDTAG_LEADING_PATTERN)) // bad formatted mdTag
		return 0;
	else
		len += boost::lexical_cast<uint32_t>(string(match[0]));
	for(searchStart = match[0].second; boost::regex_search(searchStart, mdTag.cend(), match, MDTAG_MAIN_PATTERN); searchStart = match[0].second) {
		if(match.length(1) == 1) // mismatch
			len += match.length(1);
		else
			len += match.length(1) - 1;
		len += boost::lexical_cast<uint32_t>(match.str(2));
	}
	return len;
}

uint32_t AlignmentSE::getAlnLen() const {
	uint32_t alnLen = 0;
	for(state_str::value_type s : alnPath) {
		if(bam_cigar_type(s))
			alnLen++;
	}
	return alnLen;
}

AlignmentSE& AlignmentSE::evaluate() {
	if(alnPath.empty())
		*this;
	log10P = 0;
	/* process 5' soft-clips, if any */
	for(uint64_t i = qFrom; i < alnFrom; ++i)
		log10P += (*qual)[i] / QualStr::PHRED_SCALE - (ss->clipPenalty - ss->mismatchPenalty); // use additional penalty as the difference between mismatch and clip

	for(uint64_t k = 0, i = alnFrom, j = alnStart; k < alnPath.length(); ++k) { /* k on alnPath, i on query, j on target */
		state_str::value_type s = alnPath[k];
//		fprintf(stderr, "k: %d i: %d j: %d s: %d q: %c t: %c log10P: %g\n", k, i, j, s, query->decode(i), target->decode(j), log10P);
		switch(s) {
		case BAM_CMATCH:
			if((*query)[i] & (*target)[j]) // IUPAC match (BAM_CEQUAL)
				log10P += ::log10(1 - QualStr::phredQ2P((*qual)[i]));
			else
				log10P += (*qual)[i] / QualStr::PHRED_SCALE;
			i++;
			j++;
			break;
		case BAM_CINS:
			if(k == 0 || alnPath[k - 1] != BAM_CINS) // gap open
				log10P += -ss->gapOPenalty;
			log10P += -ss->gapEPenalty;
			i++;
			break;
		case BAM_CDEL:
			if(k == 0 || alnPath[k - 1] != BAM_CDEL) // gap open
				log10P += -ss->gapOPenalty;
			log10P += -ss->gapEPenalty;
			j++;
			break;
		default:
			break;
		}
	}

	/* process 3' soft-clips, if any */
	for(uint64_t i = alnTo; i < qTo; ++i) {
		log10P += (*qual)[i] / QualStr::PHRED_SCALE - (ss->clipPenalty - ss->mismatchPenalty); // use additional penalty as the difference between mismatch and clip
	}

	return *this;
}

vector<AlignmentSE>& AlignmentSE::calcMapQ(vector<AlignmentSE>& alnList) {
	if(alnList.empty())
		return alnList;
	const size_t N = alnList.size();
//	VectorXd prior = VectorXd::Ones(N); // use a uniform prior
	VectorXd pr(N);    // alignment probability
	for(size_t i = 0; i < N; ++i)
		pr(i) = alnList[i].loglik();
	double maxV = pr.maxCoeff();
	double scale = maxV < MIN_LOGLIK_EXP ? maxV : 0;
	pr = (pr.array()-scale).exp();
	/* get postP */
//	VectorXd postP = prior.cwiseProduct(pr);
	VectorXd postP = pr / pr.sum(); // uniform prior ignored
	/* assign postP and mapQ */
	for(size_t i = 0; i < N; ++i) {
		alnList[i].postP = postP(i);
		double mapQ = QualStr::phredP2Q(1 - postP(i));
		if(::isnan(mapQ)) {
			std::cerr << "pr: " << pr.transpose() << std::endl << " postP: " << postP.transpose() << std::endl;
		}
		assert(!::isnan(mapQ));
		alnList[i].mapQ = std::min(::round(mapQ), static_cast<double> (QualStr::MAX_Q_SCORE));
	}
	return alnList;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

