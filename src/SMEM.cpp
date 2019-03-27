/*
 * SMEM.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */
#include <cassert>
#include "SMEM.h"

namespace EGriceLab {
namespace MSGseqTK {

const double SMEM::DEFAULT_MAX_EVALUE = 0.01;

SMEM& SMEM::evaluate() {
	logP = 0;
	for(int64_t i = from; i < to; ++i)
		logP += fmdidx->loglik(seq->getBase(i));
	return *this;
}

int64_t SMEM::dbDist(const SMEM& lhs, const SMEM& rhs) {
	assert(lhs.mtg == rhs.mtg);
	int64_t minD = INT64_MAX;
	for(const Loc& loc1 : lhs.locs)
		for(const Loc& loc2 : rhs.locs)
			if(isCompatitable(lhs.mtg, loc1, loc2))
				minD = std::min(minD, Loc::dist(loc1, loc2));
	return minD;
}

ostream& SMEM::write(ostream& out) const {
	/* write basic info */
	out << from << '-' << to << ':' << loglik() << ':'; /* SAstart and SAend are ignored */
	/* write Loc info */
	for(vector<GLoc>::const_iterator loc = locs.begin(); loc != locs.end(); ++loc) {
		if(loc != locs.begin())
			out << ',';
		out << *loc;
	}
	return out;
}

SMEMS SMEM::findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	SMEMS curr, prev, match;

	to = from;
	nt16_t b = seq->getBase(from);
	if(!DNAalphabet::isBasic(b)) // first base is non-basic, no-matches
		return curr;

	/* set init SMEM */
	SMEM smem(seq, mtg, fmdidx, from, from + 1, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0 = smem;
	/* forward extension */
	while(smem.to <= L) {
		if(smem.to == L) {
			curr.push_back(smem0);
			break;
		}
		else {
			smem = static_cast<const SMEM&>(smem0).fwdExt();
			if(smem.size != smem0.size) // a different [p, q, s] Bi-directional interval found
				curr.push_back(smem0);
			if(smem.size <= 0)
				break;
			// updates
			smem0 = smem;
		}
	}
	// update to
	to = smem0.to;

	/* backward extension */
	if(from == 0) {
		match = curr; // back-ext not possible
	}
	else {
		std::swap(curr, prev);
		size_t i0 = L;
		int64_t i;
		for(i = from - 1; i >= -1; --i) {
			curr.clear();
			int64_t s1 = -1;
			for(SMEMS::const_reverse_iterator smem0 = prev.rbegin(); smem0 != prev.rend(); ++smem0) { // search from the back/largest SMEM
				SMEM smem = smem0->backExt();
				if((smem.size <= 0 || i == -1) && curr.empty() && i < i0) {
					match.push_back(*smem0);
					i0 = i;
				}
				if(smem.size > 0 && s1 != smem.size) {
					curr.push_back(smem);
					s1 = smem.size;
				}
			}
			if(curr.empty())
				break;
			std::swap(curr, prev);
		}
		from = i + 1; // update from
	}
	return match;
}

SMEM& SMEM::findLocs(size_t maxNLocs) {
	locs.reserve(size);
	const size_t N = std::min(maxNLocs, static_cast<size_t>(size));
	for(size_t i = 0; i < N; ++i) {
		{
			int64_t start = fmdidx->accessSA(fwdStart + i);
			if(mtg->getStrand(start) == GLoc::FWD) // always only search loc on fwd tStrand
				locs.push_back(GLoc(start, start + length(), mtg->getLocId(start), GLoc::FWD));
		}
		{
			int64_t start = fmdidx->accessSA(revStart + i);
			if(mtg->getStrand(start) == GLoc::FWD) // always only search loc on fwd tStrand
				locs.push_back(GLoc(start, start + length(), mtg->getLocId(start), GLoc::REV));
		}
	}
	return *this;
}

int64_t SMEM::length(const SMEMS& smemChain) {
	int64_t L = 0;
	for(const SMEM& smem : smemChain)
		L += smem.length();
	return L;
}

ostream& SMEM::write(const SMEMS& smemChain, ostream& out) {
	for(SMEMS::const_iterator smem = smemChain.begin(); smem != smemChain.end(); ++smem) {
		if(smem != smemChain.begin())
			out << ';';
		out << *smem;
	}
	return out;
}

double SMEM::bestLoglik(const SMEMS& smems) {
	double minLL= inf;
	for(const SMEM& smem : smems)
		minLL = std::min(minLL, smem.loglik());
	return minLL;
}

double SMEM::bestEvalue(const SMEMS& smems) {
	double minE = inf;
	for(const SMEM& smem : smems)
		minE = std::min(minE, smem.evalue());
	return minE;
}

bool SMEM::isCompatitable(const SMEMS& smems) {
	if(smems.size() <= 1)
		return true;
	for(SMEMS::const_iterator smem = smems.begin(); smem < smems.end() - 1; ++smem)
		if(!isCompatitable(*smem, *(smem + 1)))
			return false;
	return true;
}

SMEMS SMEM::searchSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx, double maxEvalue) {
	SMEMS smemsMerged;
	double bestE = 0;
	for(int64_t from = 0, to = 0; from < seq->length(); from = bestE <= maxEvalue ? to + 1 /* good SMEMS */ : to /* bad SMEMS */) {
		SMEMS smems = SMEM::findSMEMS(seq, mtg, fmdidx, from, to);
		SMEM::filter(smems, maxEvalue); // filter bad SMEM, list may become empty
		smemsMerged.insert(smemsMerged.end(), smems.begin(), smems.end());
		bestE = SMEM::bestEvalue(smems); // update bestE
	}
	return SMEM::sort(smemsMerged);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
