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

using std::unordered_set;
const double SMEMS::DEFAULT_MAX_EVALUE = 0.01;

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
	if(!DNAalphabet::isBasic(b)) // first base is non-basic, no-matches
		return SMEM(); // return an empty SMEM

	/* init */
	SMEM smem(seq, mtg, fmdidx, from, to, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
//	assert(smem.size > 0 && smem.to <= L);
	/* forward extension */
	for(smem0 = smem; smem.to <= L && smem.size > 0; smem.fwdExt()) {
		smem0 = smem; // update smem
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
	if(!DNAalphabet::isBasic(b)) // first base is non-basic, no-matches
		return SMEM(); // return an empty SMEM

	/* init */
	SMEM smem(seq, mtg, fmdidx, from, to, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
//	assert(smem1.size > 0 && smem1.to > 0);
	/* backward extension */
	for(smem0 = smem; smem.from >= 0 && smem.size > 0; smem.backExt()) {
		smem0 = smem; // update smem
	}
	from = smem0.from;
	return smem0;
}

SMEMS SMEMS::findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	SMEMS curr, prev;

	to = from + 1;
	nt16_t b = seq->getBase(from);
	if(!DNAalphabet::isBasic(b)) // first base is non-basic, no-matches
		return curr;

	/* init SMEM */
	SMEM smem(seq, mtg, fmdidx, from, from + 1, fmdidx->getCumCount(b), fmdidx->getCumCount(DNAalphabet::complement(b)), fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b));
	SMEM smem0;
	/* forward extension */
	for(smem0 = smem; smem.to <= L; smem.fwdExt()) {
		to = smem0.to;
		if(smem.size != smem0.size) // a different BD interval found
			curr.push_back(smem0);
		if(smem.to == L) { // end found
			curr.push_back(smem);
			break;
		}
		if(smem.size <= 0)
			break;
		/* update */
		smem0 = smem;
	}

	/* backward extension of each findings */
	unordered_set<Loc> smemIdx; // a hash index for whether a SMEM has been seen (based only on from and to)
	std::swap(curr, prev);
	for(SMEM& smem : prev) {
		for(smem0 = smem; smem.from >= 0; smem.backExt()) {
			from = std::min(from, smem0.from);
			if(smem.size != smem0.size && smemIdx.count(Loc(smem0.from, smem0.to)) == 0) { // a new different BD interval found
				curr.push_back(smem0);
				smemIdx.insert(Loc(smem0.from, smem0.to));
			}
			if(smem.from == 0 && smemIdx.count(Loc(smem.from, smem.to)) == 0) { // a new begin found
				curr.push_back(smem);
				smemIdx.insert(Loc(smem.from, smem.to));
				break;
			}
			if(smem.size <= 0)
				break;
			/* update */
			smem0 = smem;
		}
	}
	// update from
	return curr;
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
