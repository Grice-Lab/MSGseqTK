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
	for(int64_t i = from; i < to; ++i) {
		DNAseq::value_type b = seq->getBase(i);
		if(strand == REV)
			b = DNAalphabet::complement(b);
		logP += std::log(fmidx->getBaseCount(b));
	}
	logP -= length() * log(fmidx->totalBases()); /* subtract denominator */
	return *this;
}

uint64_t MEM::seqDist(const MEM& mem1, const MEM& mem2) {
	assert(mem1.seq == mem2.seq);
	if(isOverlap(mem1, mem2))
		return 0;
	else
		return mem1.from < mem2.from ? mem2.from - mem1.to + 1 : mem1.from - mem2.to + 1;
}

uint64_t MEM::dbDist(const MetaGenome* mtg, const MEM& mem1, const MEM& mem2) {
	uint64_t minD = std::numeric_limits<uint64_t>::max();
	for(const Loc& loc1 : mem1.locs) {
		for(const Loc& loc2 : mem2.locs) {
			if(isCompatitable(mtg, loc1, loc2)) {
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
	out << strand << ':' << from << '-' << to << ':'; /* SAstart and SAend are ignored */
	/* write Loc info */
	for(vector<Loc>::const_iterator loc = locs.begin(); loc != locs.end(); ++loc) {
		if(loc != locs.begin())
			out << ',';
		out << *loc;
	}

	return out;
}

MEM MEM::findMEM(const PrimarySeq* seq, const FMIndex* fmidx, uint64_t from, STRAND strand) {
	uint64_t start = 0; /* 0-based SAstart */
	uint64_t end = 0;   /* 1-based SAend */
	uint64_t nextStart = start;
	uint64_t nextEnd = end;
	uint64_t to;
	/* search left-to-right */
	const DNAseq& ds = strand == FWD ? seq->getSeq() : seq->getSeq().revcom();

	for(to = from; to < ds.length(); ++to, start = nextStart, end = nextEnd) {
		sauchar_t b = ds[to];
		if(b == DNAalphabet::GAP_BASE) /* null gap */
			break;
		if(start == 0) {
			nextStart = fmidx->getCumCount(b);
			nextEnd = fmidx->getCumCount(b + 1);
		}
		else {
			nextStart = fmidx->LF(b, start - 1);
			nextEnd = fmidx->LF(b, end - 1);
		}

		if(nextStart != 0 && nextStart >= nextEnd)
			break;
	}

	/* return MEM with basic info */
	return MEM(seq, fmidx, strand, from, to, start, end);
}

MEM& MEM::findLocs(size_t maxNLocs) {
	locs.reserve(SAend - SAstart);
	for(uint64_t i = SAstart; i < std::min<size_t>(SAstart + maxNLocs, SAend); ++i) {
		uint64_t start = fmidx->accessSA(i);
		locs.push_back(fmidx->reverseLoc(Loc(start, start + length())));
//		locs.push_back(Loc(start, start + length()));
	}
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */


