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

using namespace EGriceLab::Math;
using std::unordered_set;
using std::endl;
using std::cerr;

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
	assert(from < seq->length());
	nt16_t b = seq->getBase(from);
	saidx_t p = fmdidx->getCumCount(b);
	saidx_t q = fmdidx->getCumCount(DNAalphabet::complement(b));
	saidx_t s = fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b);

	int64_t fwdStart;
	int64_t revStart;
	int64_t size;
	saidx_t to;
	for(fwdStart = p, revStart = q, size = s, to = from; s > 0 && to < seq->length() && !DNAalphabet::isGap(b); fmdidx->fwdExt(p, q, s, b)) {
		fwdStart = p;
		revStart = q;
		size = s;
		b = seq->getBase(++to);
	}
	/* return MEM with basic info */
	return MEM(seq, mtg, fmdidx, from, to, fwdStart, revStart, size);
}

MEM MEM::findMEMrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx, int64_t to) {
	assert(to > 0);
	to = std::min<int64_t>(to, seq->length());
	nt16_t b = seq->getBase(to - 1);
	saidx_t p = fmdidx->getCumCount(b);
	saidx_t q = fmdidx->getCumCount(DNAalphabet::complement(b));
	saidx_t s = fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b);

	int64_t fwdStart;
	int64_t revStart;
	int64_t size;
	saidx_t from;
	for(fwdStart = p, revStart = q, size = s, from = to; s > 0 && from > 0 && !DNAalphabet::isGap(b); fmdidx->backExt(p, q, s, b)) {
		fwdStart = p;
		revStart = q;
		size = s;
		b = seq->getBase(--from - 1);
	}
	/* return MEM with basic info */
	return MEM(seq, mtg, fmdidx, from, to, fwdStart, revStart, size);
}

MEM& MEM::findLocs(size_t maxNLocs) {
	locs.reserve(size);
	const size_t N = std::min(maxNLocs, static_cast<size_t>(size));
	for(size_t i = 0; i < N; ++i) {
		saidx_t fStart = fmdidx->accessSA(fwdStart + i);
		if(mtg->getStrand(fStart) == GLoc::FWD) // always only search loc on fwd tStrand
			locs.push_back(GLoc(fStart, fStart + length(), mtg->getLocId(fStart), GLoc::FWD));
		saidx_t rStart = fmdidx->accessSA(revStart + i);
//		std::cerr << "rStart: " << rStart << " tid: " << mtg->getLocId(rStart) << std::endl;
		if(mtg->getStrand(rStart) == GLoc::FWD) // always only search loc on fwd tStrand
			locs.push_back(GLoc(rStart, rStart + length(), mtg->getLocId(rStart), GLoc::REV));
	}
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */


