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
#include "AlignmentSE.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::cerr;
using std::endl;

const string AlignmentSE::STATES = BAM_CIGAR_STR;
const double AlignmentSE::SeedPair::DEFAULT_INDEL_RATE = 0.02;
const double AlignmentSE::DEFAULT_SEEDMATCH_EPSILON = 30;
const boost::regex AlignmentSE::MDTAG_LEADING_PATTERN = boost::regex("^\\d+");
const boost::regex AlignmentSE::MDTAG_MAIN_PATTERN = boost::regex("([A-Z]|\\^[A-Z]+)(\\d+)");

void AlignmentSE::initScores() {
//	assert(isInitiated());
	M.row(0).setZero();
	M.col(0).setZero();
	D.row(0).setConstant(infV);
	D.col(0).setConstant(infV);
	I.row(0).setConstant(infV);
	I.col(0).setConstant(infV);
}

int AlignmentSE::calculateScores(const SeedMatch& seeds) {
	assert(!seeds.empty());
	assert(seeds.isCompatitable());
	/* DP at 5' */
	calculateScores(0, seeds.front().from, 0, seeds.front().start); /* could be empty loop */

	for(SeedMatch::const_iterator seed = seeds.begin(); seed < seeds.end(); ++seed) {
		/* DP within this seed */
		calculateScores(seed->from, seed->to, seed->start, seed->end);
		/* DP in between this seed and next seed, if any */
		if(seed != seeds.end() - 1)
			calculateScores(seed->to, (seed + 1)->from, seed->end, (seed + 1)->start);
	}
	/* DP at 3' */
	calculateScores(seeds.back().to, qLen, seeds.back().end, tLen); /* could be empty loop */
	return (alnScore = M.maxCoeff(&alnTo, &alnEnd)); // determine aign 3' and score simultaneously
}

void AlignmentSE::backTrace() {
	assert(alnScore > INT_MIN);
	assert(alnTo > 0 && alnEnd > 0);
	alnPath.clear();
	alnPath.reserve(qLen + tLen); /* most time enough for back-trace */
	int32_t i = alnTo;
	int32_t j = alnEnd;
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
	alnFrom = i;
	alnStart = j;
	/* reverse alnPath */
	std::reverse(alnPath.begin(), alnPath.end());
//	assert(alnPath.front() == BAM_CMATCH && alnPath.back() == BAM_CMATCH);
}

BAM::cigar_str AlignmentSE::getAlnCigar() const {
	assert(!alnPath.empty());
	BAM::cigar_str alnCigar;
	alnCigar.reserve(alnPath.length());
	/* add 5' clip, if any */
	if(alnFrom > 0)
		alnCigar.push_back(bam_cigar_gen(alnFrom, BAM_CSOFT_CLIP));

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
	/* add 3' clip, if any */
	if(alnTo < query.length())
			alnCigar.push_back(bam_cigar_gen(query.length() - alnTo, BAM_CSOFT_CLIP));
	return alnCigar;
}

string AlignmentSE::getAlnMDTag() const {
	assert(!alnPath.empty());
	string mdTag;
	mdTag.reserve(alnPath.length());
	uint32_t prev_op = -1;
	uint32_t len = 0;
	for(uint32_t k = 0, i = alnFrom, j = alnStart; k < alnPath.length(); ++k) { // i on query, j on target, k on alnPath
		uint32_t op = alnPath[k] != BAM_CMATCH ? alnPath[k] : query[i] == target[j] ? BAM_CEQUAL : BAM_CDIFF; // differentiate = and X in MD tag
		switch(op) {
		case BAM_CEQUAL:
			i++;
			j++;
			len++;
			break;
		case BAM_CDIFF:
			mdTag.append(boost::lexical_cast<string>(len));
			len = 0;
			mdTag.push_back(DNAalphabet::decode(target[j]));
			i++;
			j++;
			break;
		case BAM_CINS:
			mdTag.append(boost::lexical_cast<string>(len));
			len = 0;
			if(op != prev_op) /* open insertion */
				mdTag.push_back(ALIGN_MD_INS);
			mdTag.push_back(DNAalphabet::decode(query[i]));
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


BAM AlignmentSE::exportBAM(const string& qname, BAM::qual_str qual, int32_t tid, int32_t tShift,
		int32_t qShift, uint16_t flag, uint8_t mapQ,
		int32_t mtid, int32_t mpos, int32_t isize, uint64_t id) const {
	assert(!alnPath.empty());
	return BAM(qname, flag, tid, tShift + alnStart - 1, mapQ,
			getAlnCigar(), qLen, nt16Encode(query), qual,
			mtid, mpos, isize, id);
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
	for(SeedMatch::size_type i = 0; i < size() - 1; ++i)
		if(SeedPair::isOverlap((*this)[i], (*this)[i + 1]))
			return false;
	return true;
}

uint32_t AlignmentSE::SeedMatch::length() const {
	uint32_t len = 0;
	for(const SeedMatch::value_type& seed : *this)
		len += seed.length();
	return len;
}

double AlignmentSE::SeedMatch::maxIndelRate() const {
	double maxIndel = 0;
	for(SeedMatch::size_type i = 0; i < size() - 1; ++i) {
		double indelRate = SeedPair::indelRate((*this)[i], (*this)[i + 1]);
		if(::abs(indelRate) > ::abs(maxIndel))
			maxIndel = indelRate;
	}
	return maxIndel;
}

AlignmentSE::SeedMatchList AlignmentSE::getSeedMatchList(const MEMS& mems, uint32_t tShift,
		uint32_t qShift, uint32_t maxIt) {
	/* get a raw SeedMatchList to store all matches of each MEMS */
	const size_t N = mems.size();
	SeedMatchList rawList, outputList;
	rawList.resize(N);
	for(size_t i = 0; i < N; ++i) {
		const MEM& mem = mems[i];
		for(const Loc& loc : mem.locs)
			rawList[i].push_back(SeedPair(mem.from - qShift, loc.start - tShift, mem.length()));
	}

	vector<size_t> idx(N, 0); // index to keep track of next element in each of the N SeedMatch
	/* non-recursive algorithm to get SeedMatchList by randomly picking up elements from rawList */
	while(outputList.size() <= maxIt) {
		/* add current combination */
		outputList.push_back(SeedMatch()); // add a default SeedMatch
		outputList.back().reserve(N);
		for(size_t i = 0; i < N; ++i)
			outputList.back().push_back(rawList[i][idx[i]]);

		/* find the rightmost SeedMatch that has more elemtns left after current element */
		size_t next = N - 1;
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

double AlignmentSE::SeedMatch::loglik() const {
	double ll = 0;
	for(const SeedPair& seed : *this)
		ll += seed.loglik();
	return ll;
}

AlignmentSE::SeedMatchList& AlignmentSE::filterSeedMatchList(SeedMatchList& sortedList, double epsilon) {
	assert(epsilon >= 0);
	double bestLoglik = sortedList.front().loglik();
	sortedList.erase(
			std::find_if(sortedList.begin(), sortedList.end(), [&](const SeedMatch& match) { return match.loglik() > bestLoglik + epsilon; }), // first element that has too large loglik
			sortedList.end());
	return sortedList;
}

/** get decoded aligned query seq */
string AlignmentSE::getAlnQSeq() const {
	string seq;
	seq.reserve(getAlnQLen());
	DNAseq::size_type i = alnFrom; // index on query
	for(state_str::value_type c : alnPath) {
		switch(c) {
		case BAM_CMATCH:
			seq.push_back(DNAalphabet::decode(query[i++]));
			break;
		case BAM_CINS:
			seq.push_back(DNAalphabet::decode(query[i++]));
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
	DNAseq::size_type i = alnStart; // index on target
	for(state_str::value_type c : alnPath) {
		switch(c) {
		case BAM_CMATCH:
			seq.push_back(DNAalphabet::decode(target[i++]));
			break;
		case BAM_CINS:
			seq.push_back(ALIGN_GAP);
			break;
		case BAM_CDEL:
			seq.push_back(DNAalphabet::decode(target[i++]));
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

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

