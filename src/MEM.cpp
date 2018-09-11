/*
 * MEM.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#include <cmath>
#include <cassert>
#include "MEM.h"
#include "Stats.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace EGriceLab::Math;

double MEM::logP() const {
	double loglik = 0;
	uint64_t N = fmidx->totalBases();
	for(uint64_t i = from; i < to; ++i)
		loglik += ::log(fmidx->getBaseCount(seq->getSeq()[i])) - ::log(N); /* use observed base frequency */
	return loglik;
}

int64_t MEM::readDist(const MEM& mem1, const MEM& mem2) {
	assert(mem1.seq == mem2.seq);
	if(isOverlap(mem1, mem2))
		return 0;
	else
		return mem1.from < mem2.from ? mem2.from - mem1.to + 1 : mem1.from - mem2.to + 1;
}

int64_t MEM::dbDist(const MEM& mem1, const MEM& mem2) {
	int64_t minD = -1;
	for(const Loc& loc1 : mem1.locs) {
		for(const Loc& loc2 : mem2.locs) {
			int64_t d = Loc::dist(loc1, loc2);
			if(minD == -1 || d < minD)
				minD = d;
		}
	}
	return minD;
}

MEM MEM::findMEM(const PrimarySeq* seq, const FMIndex* fmidx, uint64_t from) {
	uint64_t start = 0; /* 0-based SAstart */
	uint64_t end = 0;   /* 1-based SAend */
	uint64_t nextStart = start;
	uint64_t nextEnd = end;
	uint64_t to;
	/* search read left-to-right */
	for(to = from; to < seq->length(); ++to, start = nextStart, end = nextEnd) {
		sauchar_t b = seq->getSeq()[to];
		if(b == 0) /* null gap */
			break;
		if(start == 0) {
			nextStart = fmidx->getBaseCount(b);
			nextEnd = fmidx->getBaseCount(b + 1);
		}
		else {
			nextStart = fmidx->LF(b, start - 1);
			nextEnd = fmidx->LF(b, end - 1);
		}

		if(nextStart != 0 && nextStart >= nextEnd)
			break;
	}

	/* return MEM with basic info */
	return MEM(seq, fmidx, from, to, start, end);
}

MEM& MEM::findLocs() {
	for(uint64_t i = SAstart; i < SAend; ++i) {
		uint64_t start = fmidx->accessSA(i);
		locs.push_back(Loc(start, start + length()));
	}

	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */


