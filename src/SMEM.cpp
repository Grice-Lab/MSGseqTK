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

SMEM_LIST SMEM_LIST::findFwdBackSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
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
		if(smem.size != smem0.size && smem0.to >= minLen /* from begin to here is enough */
			&& (!std::isfinite(maxEvalue) || smem0.evalue(0, smem0.to) <= maxEvalue) /* from begin to here is significant */)
			prev.push_back(smem0);
		smem0 = smem;
	}
	to = smem.to;

	if(prev.empty()) {
		from--;
		return curr;
	}

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

SMEM_LIST SMEM_LIST::findBackFwdSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
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
	/* back extension */
	while(smem.size >= minSize) {
		smem.backExt();
		if(smem.size != smem0.size && L - smem0.from >= minLen /* from here to end is enough */
				&& (!std::isfinite(maxEvalue) || smem0.evalue(smem0.from, L) <= maxEvalue) /* from here to end is significant */)
			prev.push_back(smem0);
		smem0 = smem;
	}
	from = smem.from;

	if(prev.empty()) {
		to++;
		return curr;
	}

	/* forward extension */
	std::reverse(prev.begin(), prev.end()); // sort SMEM list by their length decreasingly
	while(to <= L && !prev.empty()) {
		int64_t s0 = -1; // SMEM size in prev of each iteration should be mono decreasing
		int64_t nValid = 0;
		for(SMEM& smem : prev) {
			smem0 = smem;
			smem.fwdExt();
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
		to++;
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
		int64_t minLen, double maxEvalue) {
	const int64_t L = seq->length();
	/* init SMEMS search from the mid-point */
	int64_t from = L / 2;
	int64_t to = from + 1;
	SMEM_LIST curr = SMEM_LIST::findFwdBackSMEMS(seq, mtg, fmdidx, from, to, minLen, maxEvalue);
	while(from > 0 || to < L) {
		/* fwd search */
		int64_t from0 = to;
		if(from0 < L)
			curr += SMEM_LIST::findFwdBackSMEMS(seq, mtg, fmdidx, from0, to, minLen, maxEvalue);
		/* rev search */
		int64_t to0 = --from + 1;
		if(from >= 0)
			curr += SMEM_LIST::findBackFwdSMEMS(seq, mtg, fmdidx, from, to0, minLen, maxEvalue);
	}

	/* 2nd-pass re-seeding, if the least significant SMEM is still too large */
//	SMEM_LIST prev = curr;
//	for(const SMEM& smem : prev) {
//		if(smem.length() > L - minLen && curr.totalSize() < maxNSeed) {
//			int64_t from = (smem.from + smem.to) / 2;
//			int64_t to = from + 1;
//			curr += SMEM_LIST::findFwdBackSMEMS(seq, mtg, fmdidx, from, to, minLen, maxEvalue, smem.size + 1);
//		}
//	}

	return curr;
}

SeedList SMEM::getSeeds() const {
	SeedList seeds;
	seeds.reserve(size);
	for(size_t i = 0; i < size; ++i) {
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
	assert(seeds.size() == size);
	return seeds;
}

SeedList SMEM::getSeeds(SAmap_t& SAcached) const {
	SeedList seeds;
	seeds.reserve(size);
	for(size_t i = 0; i < size; ++i) {
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
	assert(seeds.size() == size);
	return seeds;
}

SeedList SMEM_LIST::findSeeds(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t minLen, double maxEvalue, int64_t maxNSeed) {
	SMEM_LIST smems = findAllSMEMS(seq, mtg, fmdidx, minLen, maxEvalue);
	/* sort SMEMs by size */
	std::sort(smems.begin(), smems.end(),
			[](const SMEM& lhs, const SMEM& rhs)->bool { return lhs.loglik() < rhs.loglik(); });

	/* get seeds using per-seq cache */
	SeedList allSeeds;
	int64_t N = smems.totalSize();
//	if(N > maxNSeed)
//		return allSeeds;
	allSeeds.reserve(N);
	if(smems.size() > 1) { /* redundant locate may exist */
		SMEM::SAmap_t SAcached;
		for(const SMEM& smem : smems) {
			const SeedList& seeds = smem.getSeeds(SAcached);
			allSeeds.insert(allSeeds.end(), seeds.begin(), seeds.end());
			if(allSeeds.size() >= maxNSeed)
				break;
		}
	}
	else {
		for(const SMEM& smem : smems) {
			const SeedList& seeds = smem.getSeeds();
			allSeeds.insert(allSeeds.end(), seeds.begin(), seeds.end());
			if(allSeeds.size() >= maxNSeed)
				break;
		}
	}
//	assert(allSeeds.size() == N);
	/* remove redundant seeds */
	SeedPair::removeRedundant(allSeeds);
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
