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

const double SMEM_LIST::MAX_EVALUE = inf;

SMEM& SMEM::evaluate() {
	logP = 0;
	for(int64_t i = from; i < to; ++i)
		logP += fmdidx->loglik(seq->getBase(i));
	return *this;
}

double SMEM::loglik(int64_t from, int64_t to) const {
	double ll = 0;
	if(from <= this->from && to >= this->to) { // calculating in an extended range
		ll += loglik();
		for(int64_t i = from; i < this->from; ++i)
			ll += fmdidx->loglik(seq->getBase(i));
		for(int64_t i = this->to; i < to; ++i)
			ll += fmdidx->loglik(seq->getBase(i));
	}
	else { // need full evaluation
		for(int64_t i = from; i < to; ++i)
			ll += fmdidx->loglik(seq->getBase(i));
	}
	return ll;
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
		int64_t& from, int64_t& to,
		int64_t minLen, double maxEvalue, int64_t minSize) {
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
		if(smem.size != smem0.size && smem0.to >= minLen && (maxEvalue == inf || smem0.evalue(0, smem0.to) <= maxEvalue))
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
			if(smem.size != smem0.size && smem0.size >= minSize && nValid == 0 && smem0.length() >= minLen && smem0.evalue() <= maxEvalue)
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
		int64_t minLen, double maxEvalue) {
	const int64_t L = seq->length();
	MEM_LIST mems;
	for(int64_t from = 0, to = 1; from < L; from = to + 1) {
		MEM mem = SMEM::findMEM(seq, mtg, fmdidx, from, to);
		if(mem.length() >= minLen && mem.evalue() <= maxEvalue)
			mems.push_back(mem);
	}
	return mems;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen, int64_t maxLen, double maxEvalue) {
	const int64_t L = seq->length();
	SMEM_LIST curr, prev;
	/* 1st-pass SMEMS search */
	for(int64_t from = 0, to = 1; from < L; from = to + 1)
		curr += SMEM_LIST::findAllSMEMS(seq, mtg, fmdidx, from, to, minLen, maxEvalue);

	/* 2nd-pass reseeding, if requested */
	if(maxLen > 0) {
		prev = curr; // copy 1st pass result
		for(const SMEM& smem : prev) {
			if(smem.length() > maxLen) {
				int64_t from1 = (smem.from + smem.to + 1) / 2; // ceil((from + to) / 2)
				int64_t to1 = from1 + 1;
				curr += SMEM_LIST::findAllSMEMS(seq, mtg, fmdidx, from1, to1, minLen, maxEvalue, smem.size + 1);
			}
		}
	}
	std::sort(curr.begin(), curr.end(),
			[](const SMEM& lhs, const SMEM& rhs) { return lhs.loglik() < rhs.loglik(); });
	return curr;
}

SeedList SMEM::getSeeds(int64_t maxNSeed) const {
	SeedList seeds;
	const size_t N = std::min<size_t>(maxNSeed, size);
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
		int64_t minLen, int64_t maxLen, double maxEvalue, int64_t maxNSeed, bool discardSeed) {
	const SMEM_LIST& smems = findAllSMEMS(seq, mtg, fmdidx, minLen, maxLen, maxEvalue);
	/* get seeds */
	SeedList allSeeds;
	for(const SMEM& smem : smems) {
		if(discardSeed && smem.getSize() > maxNSeed)
			continue;
		const SeedList& seeds = smem.getSeeds(maxNSeed);
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
