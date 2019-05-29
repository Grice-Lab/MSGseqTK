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

const double SMEMS::DEFAULT_MAX_EVALUE = 0.001;

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

SMEMS SMEMS::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	to = from + 1;
	SMEMS curr, prev;

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
	to = getTo(curr); // update to

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
		from = getFrom(curr); // update from
		/* sort and get unique after backward extension */
		std::sort(curr.begin(), curr.end());
		curr.erase(std::unique(curr.begin(), curr.end()), curr.end());
	}
	return curr;
}

int64_t SMEMS::getFrom(const SMEMS& smems) {
	int64_t from = INT64_MAX;
	for(const SMEM& smem : smems)
		from = std::min(from, smem.getFrom());
	return from;
}

int64_t SMEMS::getTo(const SMEMS& smems) {
	int64_t to = INT64_MIN;
	for(const SMEM& smem : smems)
		to = std::max(to, smem.getTo());
	return to;
}

SMEMS SMEMS::findSMEMSfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	const int64_t L = seq->length();
	SMEMS smems;
	for(int64_t from = 0, to = 1; from < L && to <= L; from = to + 1) {
		SMEM smem = SMEM::findSMEMfwd(seq, mtg, fmdidx, from, to);
		if(smem.evalue() <= maxEvalue)
			smems.push_back(smem);
	}
	return smems;
}

SMEMS SMEMS::findSMEMSrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	const int64_t L = seq->length();
	SMEMS smems;
	for(int64_t from = L - 1, to = L; from >= 0 && to > 0; to = from - 1) {
		SMEM smem = SMEM::findSMEMrev(seq, mtg, fmdidx, from, to);
		if(smem.evalue() <= maxEvalue)
			smems.push_back(smem);
	}
	return smems;
}

SMEMS SMEMS::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		double maxEvalue) {
	SMEMS allSmems;
	for(int64_t from = 0, to = 0; from < seq->length(); from = to + 1) {
		// get SMEM list at current position
		SMEMS smems = SMEMS::findAllSMEMS(seq, mtg, fmdidx, from, to).filter(maxEvalue);
		allSmems.insert(allSmems.end(), smems.begin(), smems.end());
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
