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

const double SMEM_LIST::MAX_EVALUE = 0.001;
const double SMEM_LIST::RESEED_FACTOR = 0.25;

SMEM& SMEM::evaluate() {
	logP = 0;
	for(int64_t i = from; i < to; ++i)
		logP += fmdidx->loglik(seq->getBase(i));
	return *this;
}

SMEM SMEM::findSMEMfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	to = from + 1;
	nt16_t b = seq->getBase(from);
	if(DNAalphabet::isAmbiguous(b)) // first base is non-basic, no-matches
		return SMEM(); // return an empty SMEM

	/* init */
	SMEM smem(seq, mtg, fmdidx, from, to, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* forward extension */
	for(smem0 = smem; smem.to <= L && smem.size > 0; smem0 = smem) {
		smem.fwdExt();
	}
	to = smem0.to;
	return smem0;
}

SMEM SMEM::findSMEMrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t to) {
	const size_t L = seq->length();
	assert(to <= L);
	from = to - 1;
	nt16_t b = seq->getBase(to - 1);
	if(DNAalphabet::isAmbiguous(b)) // first base is non-basic, no-matches
		return SMEM(); // return an empty SMEM

	/* init */
	SMEM smem(seq, mtg, fmdidx, from, to, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* backward extension */
	for(smem0 = smem; smem.from >= 0 && smem.size > 0; smem0 = smem) {
		smem.backExt(); // update smem
	}
	from = smem0.from;
	return smem0;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	to = from + 1;
	SMEM_LIST curr, prev;

	nt16_t b = seq->getBase(from);
	if(DNAalphabet::isAmbiguous(b)) // first base is non-basic, no-matches
		return curr;

	/* init SMEM */
	SMEM smem(seq, mtg, fmdidx, from, from + 1, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* forward extension */
	for(smem0 = smem; smem.size > 0; smem0 = smem) {
		smem.fwdExt();
		if(smem.size != smem0.size) // a different BD interval found
			curr.push_back(smem0);
	}
	to = curr.getTo(false);

	/* backward extension, if necessary */
	if(from > 0) {
		std::swap(curr, prev);
		for(SMEM& smem : prev) {
			for(smem0 = smem; smem.size > 0; smem0 = smem) {
				smem.backExt();
				if(smem.size != smem0.size) // a different found
					curr.push_back(smem0);
			}
		}
		from = curr.getFrom(false);
	}
	return curr;
}

int64_t SMEM_LIST::getFrom(bool isSorted) const {
	int64_t from = INT64_MAX;
	if(empty())
		return from;
	if(isSorted)
		return front().getFrom();
	for(const SMEM& smem : *this)
		from = std::min(from, smem.getFrom());
	return from;
}

int64_t SMEM_LIST::getTo(bool isSorted) const {
	int64_t to = INT64_MIN;
	if(empty())
		return to;
	if(isSorted)
		return back().getTo();
	for(const SMEM& smem : *this)
		to = std::max(to, smem.getTo());
	return to;
}

SMEM_LIST SMEM_LIST::findSMEMSfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	const int64_t L = seq->length();
	SMEM_LIST smems;
	for(int64_t from = 0, to = 0; from < L; from = to + 1) {
		SMEM smem = SMEM::findSMEMfwd(seq, mtg, fmdidx, from, to);
		if(smem.evalue() <= maxEvalue)
			smems.push_back(smem);
	}
	return smems;
}

SMEM_LIST SMEM_LIST::findSMEMSrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	const int64_t L = seq->length();
	SMEM_LIST smems;
	for(int64_t from = L - 1, to = L; from >= 0 && to > 0; to = from - 1) {
		SMEM smem = SMEM::findSMEMrev(seq, mtg, fmdidx, from, to);
		if(smem.evalue() <= maxEvalue)
			smems.push_back(smem);
	}
	return smems;
}

SMEM_LIST SMEM_LIST::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	const int64_t L = seq->length();
	SMEM_LIST allSmems;
	/* 1st-pass forward search */
	for(int64_t from = 0, to = 0; from < L; from = to + 1) {
		allSmems += SMEM_LIST::findAllSMEMS(seq, mtg, fmdidx, from, to).filter(maxEvalue).sort();
	}
	return allSmems;
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
		double maxEvalue) {
	const SMEM_LIST& smems = findAllSMEMS(seq, mtg, fmdidx, maxEvalue); // get SMEMS
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
