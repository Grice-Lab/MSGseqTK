/*
 * SMEM.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */
#include <cassert>
#include "MSGseqTKConst.h"
#include "SMEM.h"

namespace EGriceLab {
namespace MSGseqTK {

const double SMEMS::DEFAULT_MAX_EVALUE = 0.01;

SMEM& SMEM::evaluate() {
	logP = 0;
	for(int64_t i = from; i < to; ++i)
		logP += fmdidx->loglik(seq->getBase(i));
	return *this;
}

SMEM SMEM::findSMEM(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	to = from + 1;
	nt16_t b = seq->getBase(from);
	if(!DNAalphabet::isBasic(b)) // first base is non-basic, no-matches
		return SMEM(); // return an empty SMEM

	/* init */
	SMEM smem(seq, mtg, fmdidx, from, to, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	assert(smem.size > 0 && smem.to <= L);
	/* forward extension */
	for(SMEM smem1 = smem; smem.to < L; smem = smem1) {
		smem1 = static_cast<const SMEM&>(smem).fwdExt();
		if(smem1.size <= 0)
			break;
	}
	to = smem.to;
	/* backward extension */
	for(SMEM smem1 = smem; smem.from > 0; smem = smem1) {
		smem1 = static_cast<const SMEM&>(smem).backExt();
		if(smem1.size <= 0)
			break;
	}
	from = smem.from;
	return smem;
}

SMEMS SMEMS::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	SMEMS curr, prev, match;

	to = from + 1;
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

SMEMS SMEMS::findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	SMEMS smems;
	SMEM smem;
	for(int64_t from = 0, to = 1; from < seq->length(); from = smem.evalue() <= maxEvalue ? to + 1 /* good SMEMS */ : to /* bad SMEMS */) {
		// get longest SMEM  at current position
		smem = SMEM::findSMEM(seq, mtg, fmdidx, from, to);
		if(smem.evalue() <= maxEvalue)
			smems.push_back(smem);
	}
	return smems;
}

SMEMS SMEMS::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	SMEMS allSmems;
	double bestE = inf;
	for(int64_t from = 0, to = 0; from < seq->length(); from = bestE <= maxEvalue ? to + 1 /* good SMEMS */ : to /* bad SMEMS */) {
		// get SMEM list at current position
		SMEMS smems = SMEMS::findAllSMEMS(seq, mtg, fmdidx, from, to);
		smems.filter(maxEvalue);
		allSmems.insert(allSmems.end(), smems.begin(), smems.end());
		bestE = std::min(bestE, smems.bestEvalue()); // update bestE
	}
	return allSmems;
}

SeedList SMEM::getSeeds() const {
	SeedList seeds;
	const size_t N = std::min<size_t>(MAX_NSEEDS, size);
	seeds.reserve(N);
	for(size_t i = 0; i < N; ++i) {
		{
			int64_t start = fmdidx->accessSA(fwdStart + i);
			int64_t tid = mtg->getLocId(start);
			if(mtg->getStrand(tid, start) == GLoc::FWD) // always only search loc on fwd tStrand
				seeds.push_back(SeedPair(from, mtg->getLoc(tid, start), length(), tid, GLoc::FWD, loglik()));
		}
		{
			int64_t start = fmdidx->accessSA(revStart + i);
			int64_t tid = mtg->getLocId(start);
			if(mtg->getStrand(tid, start) == GLoc::FWD) // always only search loc on fwd tStrand
				seeds.push_back(SeedPair(seq->length() - to, mtg->getLoc(tid, start), length(), tid, GLoc::REV, loglik()));
		}
	}
	return seeds;
}

SeedList SMEMS::findSeeds(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	const SMEMS& smems = findAllSMEMS(seq, mtg, fmdidx, maxEvalue); // get SMEMS
	/* get seeds */
	SeedList allSeeds;
	for(const SMEM& smem : smems) {
		const SeedList& seeds = smem.getSeeds();
		allSeeds.insert(allSeeds.end(), seeds.begin(), seeds.end());
	}
	/* sort seeds in lexical order */
	std::sort(allSeeds.begin(), allSeeds.end());
	return allSeeds;
}

double SMEMS::loglik() const {
	double ll = 0;
	for(const SMEM& smem : *this)
		ll += smem.loglik();
	return ll;
}

double SMEMS::bestLoglik() const {
	double minLoglik = inf;
	for(const SMEM& smem : *this)
		minLoglik = std::min(minLoglik, smem.loglik());
	return minLoglik;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
