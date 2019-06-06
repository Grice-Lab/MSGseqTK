/*
 * SMEM.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */
#include <cassert>
#include <unordered_set>
#include "MSGseqTKConst.h"
#include "SMEM.h"

namespace EGriceLab {
namespace MSGseqTK {

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
	nt16_t b = seq->getBase(from);

	/* init */
	SMEM smem(seq, mtg, fmdidx, from, from + 1, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* forward extension */
	for(smem0 = smem; smem.size > 0; smem.fwdExt())
		smem0 = smem;
	to = smem0.to;

	/* backward extension */
	for(smem0 = smem; smem.size > 0; smem.backExt())
		smem0 = smem;
	from = smem0.from;

	return smem0;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	SMEM_LIST curr, prev, match;
	nt16_t b = seq->getBase(from);

	/* init SMEM */
	SMEM smem(seq, mtg, fmdidx, from, from + 1, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* forward extension */
	for(smem0 = smem; smem.size > 0; smem0 = smem) {
		smem.fwdExt();
		if(smem.size != smem0.size) // a different BD interval found
			curr.push_back(smem0);
	}
	to = smem.to - 1;
	if(from == 0)
		return curr;

	/* backward extension */
	std::swap(curr, prev);
	std::reverse(prev.begin(), prev.end()); // keep larger SMEM near the front
	while(from >= 0 && !prev.empty()) {
		curr.clear();
		int64_t s0 = -1;
		for(const SMEM& smem : prev) {
			SMEM smem1 = smem.backExt();
			if(smem1.size == 0 && curr.empty())
				match.push_back(smem);
			if(smem1.size > 0 && smem1.size != s0) {
				s0 = smem1.size; // size already seen for this round of backExt
				curr.push_back(smem1);
			}
		}
		from--;
		std::swap(curr, prev);
	}
	return match;
}

SMEM_LIST SMEM_LIST::findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen) {
	const int64_t L = seq->length();
	SMEM_LIST smems;
	for(int64_t from = 0, to = 1; from < L; from = to + 1) {
		SMEM smem = SMEM::findSMEM(seq, mtg, fmdidx, from, to);
		if(smem.length() >= minLen)
			smems.push_back(smem);
	}
	return smems;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen) {
	const int64_t L = seq->length();
	SMEM_LIST allSmems;
	/* forward SMEMS at 5' */
	for(int64_t from = 0, to = 1; from < L; from = to + 1) {
		allSmems += SMEM_LIST::findAllSMEMS(seq, mtg, fmdidx, from, to).filter(minLen).sort();
	}
	return allSmems.sort();
}

SeedList SMEM::getSeeds() const {
	SeedList seeds;
	const size_t N = std::min<size_t>(MAX_NSEEDS, size);
	seeds.reserve(N);
	for(size_t i = 0; i < N; ++i) {
		{
			int64_t bdStart = fmdidx->accessSA(fwdStart + i);
			int64_t tid = mtg->getTid(bdStart);
			int64_t start = mtg->getLoc(tid, bdStart);
			assert(tid == mtg->getTid(bdStart + length()));
			if(mtg->getStrand(tid, bdStart) == GLoc::FWD) { // always only search loc on fwd tStrand
				seeds.push_back(SeedPair(from, start, length(), tid, GLoc::FWD, loglik()));
			}
		}
		{
			int64_t bdStart = fmdidx->accessSA(revStart + i);
			int64_t tid = mtg->getTid(bdStart);
			int64_t start = mtg->getLoc(tid, bdStart);
			if(mtg->getStrand(tid, bdStart) == GLoc::FWD) { // always only search loc on fwd tStrand
				assert(tid == mtg->getTid(bdStart + length()));
				seeds.push_back(SeedPair(seq->length() - to, start, length(), tid, GLoc::REV, loglik()));
			}
		}
	}
	return seeds;
}

SeedList SMEM_LIST::findSeeds(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen) {
	const SMEM_LIST& smems = findAllSMEMS(seq, mtg, fmdidx, minLen);
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

double SMEM_LIST::loglik() const {
	double ll = 0;
	for(const SMEM& smem : *this)
		ll += smem.loglik();
	return ll;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
