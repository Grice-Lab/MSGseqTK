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
#include "MEM.h"
#include "Stats.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace EGriceLab::Math;
using std::unordered_set;

double MEM::loglik() const {
	double loglik = 0;
	uint64_t N = fmidx->totalBases();
	DNAseq ds = seq->getSeq();
	if(strand == REV)
		ds.revcom();
	for(uint64_t i = from; i < to; ++i)
		loglik += ::log(fmidx->getBaseCount(ds[i])) - ::log(N); /* use observed base frequency */
	return loglik;
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
	DNAseq ds = seq->getSeq();
	if(strand == REV)
		ds.revcom();

	for(to = from; to < ds.length(); ++to, start = nextStart, end = nextEnd) {
		sauchar_t b = ds[to];
		if(b == DNAalphabet::N) /* null gap */
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

MEM& MEM::findLocs() {
	for(uint64_t i = SAstart; i < SAend; ++i) {
		uint64_t start = fmidx->accessSA(i);
		locs.push_back(Loc(start, start + length()));
	}
	return *this;
}

void MEM::filterLocsByIndel(const MetaGenome* mtg, MEM& mem1, MEM& mem2, double maxIndelRate) {
	assert(mem1.seq == mem2.seq);
	unordered_set<size_t> acceptedLocs1;
	unordered_set<size_t> acceptedLocs2;
	for(size_t i = 0; i < mem1.locs.size(); ++i)
		for(size_t j = 0; j < mem2.locs.size(); ++j)
			if(isCompatitable(mtg, mem1.locs[i], mem2.locs[j])) {
				int64_t d = Loc::dist(mem1.locs[i], mem2.locs[j]);
				if(::abs(static_cast<double>(d) / mem1.seq->length() <= maxIndelRate)) {
					acceptedLocs1.insert(i);
					acceptedLocs2.insert(j);
				}
			}

	/* filter mem1 locs */
	for(size_t i = mem1.locs.size(); i > 0; --i)
		if(acceptedLocs1.count(i - 1) == 0) /* not acceptible */
			mem1.locs.erase(mem1.locs.begin() + i - 1);

	/* filter mem2 locs */
	for(size_t j = mem2.locs.size(); j > 0; --j)
		if(acceptedLocs2.count(j - 1) == 0) /* not acceptible */
			mem2.locs.erase(mem2.locs.begin() + j - 1);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */


