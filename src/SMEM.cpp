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

MEM SMEM::findMEM(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	nt16_t b = seq->getBase(from);

	/* init */
	MEM mem(seq, mtg, fmdidx, from, from + 1, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	MEM mem0;
	/* forward extension */
	for(mem0 = mem; mem.size > 0; mem.fwdExt())
		mem0 = mem;
	to = mem0.to;

	/* backward extension */
	for(mem0 = mem; mem.size > 0; mem.backExt())
		mem0 = mem;
	from = mem0.from;
	return mem0;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to, int64_t minLen, int64_t minSize) {
	const size_t L = seq->length();
	assert(from < L);
	SMEM_LIST curr, prev;
	nt16_t b = seq->getBase(from);
	to = from + 1;
	if(DNAalphabet::isAmbiguous(b))
		return curr;

	/* init SMEM */
	SMEM smem(seq, mtg, fmdidx, from, to, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* forward extension */
	for(smem0 = smem; smem.size >= minSize; smem0 = smem) {
		smem.fwdExt();
		if(smem.size != smem0.size && smem0.to >= minLen)
			prev.push_back(smem0);
	}
	to = smem.to - 1;

	/* backward extension */
	std::reverse(prev.begin(), prev.end()); // sort SMEM list by their length decreasingly
	while(from >= 0 && !prev.empty()) {
		int64_t s0 = -1; // SMEM size in prev of each iteration should be mono decreasing
		int64_t nValid = 0;
		for(SMEM& smem : prev) {
			smem0 = smem;
			smem.backExt();
			if(smem.size != smem0.size && nValid == 0 && smem0.length() >= minLen)
				curr.push_back(smem0);
			if(smem.size >= minSize && smem.size != s0) { /* still a valid SMEM with size not seen */
				s0 = smem.size;
				nValid++;
			}
			else // mark this SMEM as not valid
				smem.size = 0;
		}
		/* remove invalid SMEM */
		prev.erase(std::remove_if(prev.begin(), prev.end(),
				[](const SMEM& smem) { return !smem.isValid(); }), prev.end());
		assert(prev.size() == nValid);
		from--;
	}
	return curr;
}

MEM_LIST SMEM_LIST::findMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen) {
	const int64_t L = seq->length();
	MEM_LIST mems;
	for(int64_t from = 0, to = 1; from < L; from = to + 1) {
		MEM mem = SMEM::findMEM(seq, mtg, fmdidx, from, to);
		if(mem.length() >= minLen)
			mems.push_back(mem);
	}
	return mems;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen, int64_t maxLen) {
	const int64_t L = seq->length();
	SMEM_LIST allSmems;
	/* 1st-pass SMEMS search */
	for(int64_t from = 0, to = 1; from < L; from = to + 1)
		allSmems += SMEM_LIST::findAllSMEMS(seq, mtg, fmdidx, from, to, minLen);

	/* 2nd-pass reseeding, if requested */
	if(maxLen > 0) {
		for(const SMEM& smem : allSmems) {
			if(smem.length() > maxLen) {
				int64_t from1 = (smem.from + smem.to) / 2 - 1;
				int64_t to1 = from1 + 1;
				allSmems += SMEM_LIST::findAllSMEMS(seq, mtg, fmdidx, from1, to1, minLen, smem.size + 1);
			}
		}
	}
	return allSmems;
}

SeedList SMEM::getSeeds() const {
	SeedList seeds;
	const size_t N = std::min<size_t>(MAX_NSEED, size);
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
		int64_t minLen, int64_t maxLen) {
	const SMEM_LIST& smems = findAllSMEMS(seq, mtg, fmdidx, minLen, maxLen);
	/* get seeds */
	SeedList allSeeds;
	for(const SMEM& smem : smems) {
		const SeedList& seeds = smem.getSeeds();
		allSeeds.insert(allSeeds.end(), seeds.begin(), seeds.end());
	}
	/* sort and get unique seeds */
	std::sort(allSeeds.begin(), allSeeds.end());
	allSeeds.erase(std::unique(allSeeds.begin(), allSeeds.end()), allSeeds.end());
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
