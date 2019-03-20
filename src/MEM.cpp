/*
 * MEM.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#include <cctype>
#include <cmath>
#include <cassert>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include "MEM.h"
#include "Stats.h"

namespace EGriceLab {
namespace MSGseqTK {

MEM& MEM::evaluate() {
	if(empty())
		return *this;
	logP = 0;
	int64_t N = 0;
	for(int64_t i = from; i < to; ++i) {
		DNAseq::value_type b = seq->getBase(i);
		if(fmdidx->getBaseCount(b) > 0) { // ignore non-existing bases
			logP += std::log(fmdidx->getBaseCount(b));
			N++;
		}
	}
	logP -= N * log(fmdidx->length()); /* subtract denominator */
	return *this;
}

int64_t MEM::dbDist(const MEM& lhs, const MEM& rhs) {
	assert(lhs.mtg == rhs.mtg);
	int64_t minD = std::numeric_limits<int64_t>::max();
	for(const Loc& loc1 : lhs.locs) {
		for(const Loc& loc2 : rhs.locs) {
			if(isCompatitable(lhs.mtg, loc1, loc2)) {
				int64_t d = Loc::dist(loc1, loc2);
				if(d < minD)
					minD = d;
			}
		}
	}
	return minD;
}

ostream& MEM::write(ostream& out) const {
	/* write basic info */
	out << from << '-' << to << ':'; /* SAstart and SAend are ignored */
	/* write Loc info */
	for(vector<GLoc>::const_iterator loc = locs.begin(); loc != locs.end(); ++loc) {
		if(loc != locs.begin())
			out << ',';
		out << *loc;
	}
	return out;
}

MEM MEM::findMEMfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx, int64_t from) {
	const size_t L = seq->length();
	assert(from < L);
	int64_t to = from;
	nt16_t b = seq->getBase(to);
	int64_t p = fmdidx->getCumCount(b);
	int64_t q = fmdidx->getCumCount(DNAalphabet::complement(b));
	int64_t s = fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b);
	for(to = from + 1; to < L && fmdidx->fwdExt(p, q, s, seq->getBase(to)); ++to)
		continue;
	return MEM(seq, mtg, fmdidx, from, to, p, q, s);
}

MEM MEM::findMEMrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx, int64_t to) {
	const size_t L = seq->length();
	to = std::min<int64_t>(to, L);
	assert(0 < to && to <= L);
	int64_t from = to;
	nt16_t b = seq->getBase(from - 1);
	int64_t p = fmdidx->getCumCount(b);
	int64_t q = fmdidx->getCumCount(DNAalphabet::complement(b));
	int64_t s = fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b);
	for(from = to - 1; from > 0 && fmdidx->backExt(p, q, s, seq->getBase(from - 1)); --from)
		continue;
	return MEM(seq, mtg, fmdidx, from, to, p, q, s);
}

MEM& MEM::findLocs(size_t maxNLocs) {
	locs.reserve(size);
	const size_t N = std::min(maxNLocs, static_cast<size_t>(size));
	for(size_t i = 0; i < N; ++i) {
		{
			saidx_t start = fmdidx->accessSA(fwdStart + i);
			if(mtg->getStrand(start) == GLoc::FWD) // always only search loc on fwd tStrand
				locs.push_back(GLoc(start, start + length(), mtg->getLocId(start), GLoc::FWD));
		}
		{
			saidx_t start = fmdidx->accessSA(revStart + i);
			if(mtg->getStrand(start) == GLoc::FWD) // always only search loc on fwd tStrand
				locs.push_back(GLoc(start, start + length(), mtg->getLocId(start), GLoc::REV));
		}
	}
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */


