/*
 * SMEM.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */

#include <cmath>
#include <cassert>
#include <unordered_set>
#include "MSGseqTKConst.h"
#include "SMEM.h"

namespace EGriceLab {
namespace MSGseqTK {

const double SMEM::MAX_EVALUE = 0.01;
static const int64_t MIN_SA_CACHE_SIZE = 2; /* minimum SMEMS size to trigger using the cashe findSeeds */

SMEM& SMEM::evaluate() {
	logP = 0;
	for(int64_t i = from; i < to; ++i)
		logP += fmdidx->loglik(seq->getBase(i));
	return *this;
}

double SMEM::loglik(int64_t from, int64_t to) const {
	double ll = 0;
	if(from <= this->from && this->to <= to) { // calculating in an extended range
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
	MEM mem0(mem);
	/* forward extension */
	while(mem.fwdExt().size > 0)
			mem0 = mem;
	to = mem0.to;

	/* backward extension */
	mem0 = mem;
	while(mem.backExt().size > 0)
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
	SMEM smem0 = smem;
	/* forward extension */
	while(smem.size >= minSize) {
		smem.fwdExt();
		if(smem.size != smem0.size && smem0.to >= minLen
				&& (!std::isfinite(maxEvalue) || smem0.evalue(0, smem0.to) <= maxEvalue))
			prev.push_back(smem0);
		smem0 = smem;
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
//			std::cerr << "backexting smem0: " << smem0 << " smem: " << smem << std::endl;
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
	seeds.reserve(std::min(maxNSeed, size));
	for(size_t i = 0; i < size && seeds.size() < maxNSeed; ++i) {
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
			assert(tid == mtg->getTid(bdStart + length()));
			if(mtg->getStrand(tid, bdStart) == GLoc::FWD) { // always only search loc on fwd tStrand
				seeds.push_back(SeedPair(seq->length() - to, start, length(), tid, GLoc::REV, loglik()));
			}
		}
	}
	return seeds;
}

SeedList SMEM::getSeeds(SAmap_t& SAcached, int64_t maxNSeed) const {
	SeedList seeds;
	seeds.reserve(std::min(maxNSeed, size));
	for(size_t i = 0; i < size && seeds.size() < maxNSeed; ++i) {
		{
			int64_t bdStart = SAcached.count(fwdStart + i) > 0 ? SAcached.at(fwdStart + i) : (SAcached[fwdStart + i] = fmdidx->accessSA(fwdStart + i));
			int64_t tid = mtg->getTid(bdStart);
			int64_t start = mtg->getLoc(tid, bdStart);
			assert(tid == mtg->getTid(bdStart + length()));
			if(mtg->getStrand(tid, bdStart) == GLoc::FWD) { // always only search loc on fwd tStrand
				seeds.push_back(SeedPair(from, start, length(), tid, GLoc::FWD, loglik()));
			}
		}
		{
			int64_t bdStart = SAcached.count(revStart + i) > 0 ? SAcached.at(revStart + i) : (SAcached[revStart + i] = fmdidx->accessSA(revStart + i));
			int64_t tid = mtg->getTid(bdStart);
			int64_t start = mtg->getLoc(tid, bdStart);
			assert(tid == mtg->getTid(bdStart + length()));
			if(mtg->getStrand(tid, bdStart) == GLoc::FWD) { // always only search loc on fwd tStrand
				seeds.push_back(SeedPair(seq->length() - to, start, length(), tid, GLoc::REV, loglik()));
			}
		}
	}
	return seeds;
}

SeedList SMEM_LIST::findSeeds(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen, int64_t maxLen, double maxEvalue, int64_t maxNSeed) {
	const SMEM_LIST& smems = findAllSMEMS(seq, mtg, fmdidx, minLen, maxLen, maxEvalue);
	/* get seeds using per-seq cache */
	SeedList allSeeds;
	allSeeds.reserve(maxNSeed * smems.size());
	if(smems.size() < MIN_SA_CACHE_SIZE) {
		for(const SMEM& smem : smems) {
			const SeedList& seeds = smem.getSeeds(maxNSeed);
			allSeeds.insert(allSeeds.end(), seeds.begin(), seeds.end());
		}
	}
	else {
		SMEM::SAmap_t SAcached;
		for(const SMEM& smem : smems) {
			const SeedList& seeds = smem.getSeeds(SAcached, maxNSeed);
			allSeeds.insert(allSeeds.end(), seeds.begin(), seeds.end());
		}
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
