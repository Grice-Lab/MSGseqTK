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

const double MEMS::DEFAULT_MAX_EVALUE = 0.01;

uint64_t MEMS::length() const {
	uint64_t L = 0;
	for(const MEM& mem : *this)
		L += mem.length();
	return L;
}

MEMS MEMS::sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx, RNG& rng, double maxEvalue,
		uint64_t from, uint64_t to, MEM::STRAND strand) {
	MEMS mems;
	for(uint64_t i = from; i < std::min(to, seq->length());) {
		MEM mem = MEM::findMEM(seq, fmidx, i, to, strand).evaluate();

		/* calculate MEM evalue */
		double eval = mem.evalue();
		/* accept by chance */
		if(eval <= maxEvalue && eval <= mem_dist(rng)) {
			mems.push_back(mem);
			i = mem.to + 1;
		}
		else
			i = eval <= mem_dist(rng) ? mem.to + 1 : i + 1;
	}
	return mems;
}

MEMS MEMS::sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx, RNG& rng, double maxEvalue,
		uint64_t from, uint64_t to, int strand) {
	assert(strand != 0);
	if((strand & MEM::FWD) && !(strand & MEM::REV)) /* fwd only */
		return sampleMEMS(seq, fmidx, rng, maxEvalue, from, to, MEM::FWD);
	else if((strand & MEM::FWD) && !(strand & MEM::REV)) /* rev only */
		return sampleMEMS(seq, fmidx, rng, maxEvalue, from, to, MEM::REV);
	else {
		MEMS fwdMems = sampleMEMS(seq, fmidx, rng, maxEvalue, from, to, MEM::FWD);
		MEMS revMems = sampleMEMS(seq, fmidx, rng, maxEvalue, from, to, MEM::REV);
		return fwdMems.loglik() < revMems.loglik() ? fwdMems : revMems; /* return the most significant result */
	}
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

MEMS_PE MEMS_PE::sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx, RNG& rng, double maxEvalue,
		uint64_t fwdFrom, uint64_t fwdTo, uint64_t revFrom, uint64_t revTo, int strand) {
	assert(strand != 0);
	if((strand & MEM::FWD) && !(strand & MEM::REV)) /* fwd only */
		return MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, maxEvalue, fwdFrom, fwdTo, MEM::FWD),
				MEMS::sampleMEMS(revSeq, fmidx, rng, maxEvalue, revFrom, revTo, MEM::REV)); // use FWD-REV orientation
	else if(!(strand & MEM::FWD) && (strand & MEM::REV)) /* rev only */
		return MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, maxEvalue, fwdFrom, fwdTo, MEM::REV),
				MEMS::sampleMEMS(revSeq, fmidx, rng, maxEvalue, revFrom, revTo, MEM::FWD)); // use REV-FWD orientation
	else {
		MEMS_PE senseMemsPE = MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, maxEvalue, fwdFrom, fwdTo, MEM::FWD),
				MEMS::sampleMEMS(revSeq, fmidx, rng, maxEvalue, revFrom, revTo, MEM::REV)); // use FWD-REV orientation
		MEMS_PE antiMemsPE = MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, maxEvalue, fwdFrom, fwdTo, MEM::REV),
				MEMS::sampleMEMS(revSeq, fmidx, rng, maxEvalue, revFrom, revTo, MEM::FWD)); // use REV-FWD orientation
		return senseMemsPE.loglik() < antiMemsPE.loglik() ? senseMemsPE : antiMemsPE; /* return the most significant result */
	}
}

} /* namespace UCSC */
} /* namespace EGriceLab */
