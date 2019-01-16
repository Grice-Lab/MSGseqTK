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

MEMS MEMS::sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
		RNG& rng, int strand) {
	assert(strand != 0);
	if((strand & MEM::FWD) && !(strand & MEM::REV)) /* fwd only */
		return sampleMEMS(seq, fmidx, rng, MEM::FWD);
	else if((strand & MEM::FWD) && !(strand & MEM::REV)) /* rev only */
		return sampleMEMS(seq, fmidx, rng, MEM::REV);
	else {
		MEMS fwdMems = sampleMEMS(seq, fmidx, rng, MEM::FWD);
		MEMS revMems = sampleMEMS(seq, fmidx, rng, MEM::REV);
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

MEMS_PE MEMS_PE::sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx,
		RNG& rng, int strand) {
	assert(strand != 0);
	if((strand & MEM::FWD) && !(strand & MEM::REV)) /* fwd only */
		return MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, MEM::FWD),
				MEMS::sampleMEMS(revSeq, fmidx, rng, MEM::REV)); // use FWD-REV orientation
	else if(!(strand & MEM::FWD) && (strand & MEM::REV)) /* rev only */
		return MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, MEM::REV),
				MEMS::sampleMEMS(revSeq, fmidx, rng, MEM::FWD)); // use REV-FWD orientation
	else {
		MEMS_PE senseMemsPE = MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, MEM::FWD),
				MEMS::sampleMEMS(revSeq, fmidx, rng, MEM::REV)); // use FWD-REV orientation
		MEMS_PE antiMemsPE = MEMS_PE(MEMS::sampleMEMS(fwdSeq, fmidx, rng, MEM::REV),
				MEMS::sampleMEMS(revSeq, fmidx, rng, MEM::FWD)); // use REV-FWD orientation
		return senseMemsPE.loglik() < antiMemsPE.loglik() ? senseMemsPE : antiMemsPE; /* return the most significant result */
	}
}

} /* namespace UCSC */
} /* namespace EGriceLab */
