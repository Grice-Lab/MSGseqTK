/*
 * MEMS.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: zhengqi
 */

#include <unordered_set>
#include <stack>
#include <cmath>
#include <utility>
#include "MEMS.h"
#include "MSGseqTKConst.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace EGriceLab {
namespace MSGseqTK {

using std::unordered_set;
using std::stack;

uint64_t MEMS::length() const {
	uint64_t L = 0;
	for(const MEM& mem : *this)
		L += mem.length();
	return L;
}

MEMS MEMS::sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
		RNG& rng, MEM::STRAND strand) {
	MEMS mems;
	for(int64_t i = 0; i < seq->length();) {
		MEM mem = MEM::findMEM(seq, fmidx, i, strand).evaluate();

		/* calculate MEM evalue */
		double eval = mem.evalue();
		/* accept by chance */
		if(mem.evalue() <= mem_dist(rng)) {
			mems.push_back(mem);
			i = mem.to + 1;
		}
		else
			i++;
	}
	return mems;
}

MEMS MEMS::getBestMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
		RNG& rng, int strand, int nSeed) {
	assert(strand != 0);
	MEMS fwdMems, revMems; // best fwd/rev MEMS
	if(strand & MEM::FWD) { // FWD strand need search
		double bestLogP = infV;
		for(int k = 0; k < nSeed; ++k) {
			MEMS mems = sampleMEMS(seq, fmidx, rng, MEM::FWD);
			if(mems.loglik() > bestLogP) {
				fwdMems = mems;
				bestLogP = mems.loglik();
			}
		}
	}
	if(strand & MEM::REV) { // FWD strand need search
		double bestLogP = infV;
		for(int k = 0; k < nSeed; ++k) {
			MEMS mems = sampleMEMS(seq, fmidx, rng, MEM::REV);
			if(mems.loglik() > bestLogP) {
				revMems = mems;
				bestLogP = mems.loglik();
			}
		}
	}

	if(strand & MEM::FWD != 0 && strand & MEM::REV == 0) /* fwd only */
		return fwdMems;
	else if(strand & MEM::FWD == 0 && strand & MEM::REV != 0)
		return revMems;
	else
		return fwdMems.loglik() < revMems.loglik() ? fwdMems : revMems; /* return the most significant result */
}

size_t MEMS::bestMEMIndex() const {
	size_t bestIdx = -1;
	double maxLoglik = inf;
	for(MEMS::size_type i = 0; i < size(); ++i)
		if((*this)[i].loglik() < maxLoglik)
			bestIdx = i;
	return bestIdx;
}

ostream& MEMS::write(ostream& out) const {
	for(const MEM& mem : *this)
		out << mem << ';';
	out << "strand=" << getStrand() << ";loglik=" << loglik();

	return out;
}

MEMS_PE MEMS::getBestMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx,
		RNG& rng, int strand, int nSeed) {
	assert(strand != 0);
	MEMS_PE senseMemsPE, antiMemsPE; // sense and anti-sense MEMS_PE
	/* scan FWD strand */
	if(strand & MEM::FWD) {
		double bestLogP = infV;
		for(int k = 0; k < nSeed; ++k) {
			MEMS_PE memsPE = std::make_pair(sampleMEMS(fwdSeq, fmidx, rng, MEM::FWD),
					sampleMEMS(revSeq, fmidx, rng, MEM::REV)); // use FWD-REV orientation
			if(loglik(memsPE) > bestLogP) {
				senseMemsPE = memsPE;
				bestLogP = loglik(memsPE);
			}
		}
	}
	/* scan REV strand */
	if(strand & MEM::REV) {
		double bestLogP = infV;
		for(int k = 0; k < nSeed; ++k) {
			MEMS_PE memsPE = std::make_pair(sampleMEMS(fwdSeq, fmidx, rng, MEM::REV),
					sampleMEMS(revSeq, fmidx, rng, MEM::FWD)); // use REV-FWD orientation
			if(loglik(memsPE) > bestLogP) {
				antiMemsPE = memsPE;
				bestLogP = loglik(memsPE);
			}
		}
	}

	if(strand & MEM::FWD != 0 && strand & MEM::REV == 0) /* fwd only */
		return senseMemsPE;
	else if(strand & MEM::FWD == 0 && strand & MEM::REV != 0)
		return antiMemsPE;
	else
		return loglik(senseMemsPE) < loglik(antiMemsPE) ? senseMemsPE : antiMemsPE; /* return the most significant result */
}

} /* namespace UCSC */
} /* namespace EGriceLab */
