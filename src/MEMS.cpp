/*
 * MEMS.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: zhengqi
 */

#include <unordered_set>
#include <stack>
#include <cmath>
#include "MEMS.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::unordered_set;
using std::stack;

void MEMS::evaluate() {
	logP = 0;
	for(const MEM& mem : *this)
		logP += mem.loglik();
}

uint64_t MEMS::length() const {
	uint64_t L = 0;
	for(const MEM& mem : *this)
		L += mem.length();
	return L;
}

MEMS MEMS::sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
		RNG& rng, int strand, bool keepLoc) {
	assert(strand != 0);
	MEMS fwdMems, revMems;
	/* scan fwd strand */
	if(strand & MEM::FWD != 0) {
		for(int64_t i = 0; i < seq->length();) {
			MEM& mem = MEM::findMEM(seq, fmidx, i, MEM::FWD).evaluate();

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				fwdMems.push_back(!keepLoc ? mem : mem.findLocs());
				i = mem.to + 1;
			}
			else
				i++;
		}
		fwdMems.evaluate();
	}
	/* scal rev strand */
	if(strand & MEM::REV != 0) {
		for(int64_t i = 0; i < seq->length();) {
			MEM& mem = MEM::findMEM(seq, fmidx, i, MEM::REV).evaluate();
			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				revMems.push_back(!keepLoc ? mem : mem.findLocs());
				i = mem.to + 1;
			}
			else
				i++;
		}
		revMems.evaluate();
	}

	if(strand & MEM::FWD != 0 && strand & MEM::REV == 0) /* fwd only */
		return fwdMems;
	else if(strand & MEM::FWD == 0 && strand & MEM::REV != 0) /* rev only */
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

MEMS_PE MEMS::sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx,
		RNG& rng, int strand, bool keepLoc) {
	assert(strand != 0);
	MEMS_PE sense_mems_pe, revcom_mems_pe;
	/* scan fwd strand */
	if(strand & MEM::FWD != 0) {
		for(int64_t i = 0; i < fwdSeq->length();) {
			MEM& mem = MEM::findMEM(fwdSeq, fmidx, i, MEM::FWD).evaluate();

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				sense_mems_pe.first.push_back(!keepLoc ? mem : mem.findLocs());
				i = mem.to + 1;
			}
			else
				i++;
		}
		for(int64_t i = 0; i < revSeq->length();) {
			MEM& mem = MEM::findMEM(revSeq, fmidx, i, MEM::REV).evaluate();

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				sense_mems_pe.second.push_back(!keepLoc ? mem : mem.findLocs());
				i = mem.to + 1;
			}
			else
				i++;
		}
		sense_mems_pe.first.evaluate();
		sense_mems_pe.second.evaluate();
	}
	/* scal rev strand */
	if(strand & MEM::REV != 0) {
		for(int64_t i = 0; i < fwdSeq->length();) {
			MEM& mem = MEM::findMEM(fwdSeq, fmidx, i, MEM::REV).evaluate();

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				revcom_mems_pe.first.push_back(!keepLoc ? mem : mem.findLocs());
				i = mem.to + 1;
			}
			else
				i++;
		}
		for(int64_t i = 0; i < revSeq->length();) {
			MEM& mem = MEM::findMEM(revSeq, fmidx, i, MEM::FWD).evaluate();
			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				revcom_mems_pe.second.push_back(!keepLoc ? mem : mem.findLocs());
				i = mem.to + 1;
			}
			else
				i++;
		}
		sense_mems_pe.first.evaluate();
		sense_mems_pe.second.evaluate();
	}

	if(strand & MEM::FWD != 0 && strand & MEM::REV == 0) /* fwd only */
		return sense_mems_pe;
	else if(strand & MEM::FWD == 0 && strand & MEM::REV != 0)
		return revcom_mems_pe;
	else
		return loglik(sense_mems_pe) < loglik(revcom_mems_pe) ? sense_mems_pe : revcom_mems_pe; /* return the most significant result */
}

} /* namespace UCSC */
} /* namespace EGriceLab */
