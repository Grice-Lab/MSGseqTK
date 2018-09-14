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

double MEMS::loglik() const {
	double ll = 0;
	for(const MEM& mem : *this)
		ll += mem.loglik();
	return ll;
}

MEMS MEMS::sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx, const MetaGenome* mtg,
		RNG& rng, int strand, double maxIndelRate) {
	assert(strand != 0);
	MEMS fwdMems, revMems;
	double fwdLoglik = 0;
	double revLoglik = 0;
	/* scan fwd strand */
	if(strand & MEM::FWD != 0) {
		for(int64_t i = 0; i < seq->length();) {
			MEM mem = MEM::findMEM(seq, fmidx, i, MEM::FWD);

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = false;
			if(eval <= mem_dist(rng)) {
				mem.findLocs();
				acceptible = fwdMems.empty() || ::abs(MEM::rIndel(mtg, fwdMems.back(), mem)) <= maxIndelRate;
			}
			if(acceptible) {
				fwdMems.push_back(mem);
				i = mem.to + 1;
				fwdLoglik += mem.loglik();
			}
			else
				i++;
		}
	}
	/* scal rev strand */
	if(strand & MEM::REV != 0) {
		for(int64_t i = 0; i < seq->length();) {
			MEM mem = MEM::findMEM(seq, fmidx, i, MEM::REV);
			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = false;
			if(eval <= mem_dist(rng)) {
				mem.findLocs();
				acceptible = revMems.empty() || ::abs(MEM::rIndel(mtg, revMems.back(), mem)) <= maxIndelRate;
			}
			if(acceptible) {
				revMems.push_back(mem);
				i = mem.to + 1;
				revLoglik += mem.loglik();
			}
			else
				i++;
		}
	}

	if(strand & MEM::FWD != 0 && strand & MEM::REV == 0) /* fwd only */
		return fwdMems;
	else if(strand & MEM::FWD == 0 && strand & MEM::REV != 0)
		return revMems;
	else
		return fwdLoglik < revLoglik ? fwdMems : revMems; /* return the most significant result */
}

size_t MEMS::bestMEMIndex() const {
	size_t bestIdx = -1;
	double maxLoglik = inf;
	for(MEMS::size_type i = 0; i < size(); ++i)
		if((*this)[i].loglik() < maxLoglik)
			bestIdx = i;
	return bestIdx;
}

MEMS_PE MEMS::sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx, const MetaGenome* mtg,
		RNG& rng, int strand, double maxIndelRate) {
	assert(strand != 0);
	MEMS_PE sense_mems_pe, revcom_mems_pe;
	double senseLoglik = 0;
	double revcomLoglik = 0;
	/* scan fwd strand */
	if(strand & MEM::FWD != 0) {
		for(int64_t i = 0; i < fwdSeq->length();) {
			MEM mem = MEM::findMEM(fwdSeq, fmidx, i, MEM::FWD);

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = false;
			if(eval <= mem_dist(rng)) {
				mem.findLocs();
				acceptible = sense_mems_pe.first.empty() || ::abs(MEM::rIndel(mtg, sense_mems_pe.first.back(), mem)) <= maxIndelRate;
			}
			if(acceptible) {
				sense_mems_pe.first.push_back(mem);
				i = mem.to + 1;
				senseLoglik += mem.loglik();
			}
			else
				i++;
		}
		for(int64_t i = 0; i < revSeq->length();) {
			MEM mem = MEM::findMEM(revSeq, fmidx, i, MEM::REV);

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = false;
			if(eval <= mem_dist(rng)) {
				mem.findLocs();
				acceptible = sense_mems_pe.second.empty() || ::abs(MEM::rIndel(mtg, sense_mems_pe.second.back(), mem)) <= maxIndelRate;
				sense_mems_pe.second.push_back(mem);
				i = mem.to + 1;
				senseLoglik += mem.loglik();
			}
			else
				i++;
		}
	}
	/* scal rev strand */
	if(strand & MEM::REV != 0) {
		for(int64_t i = 0; i < fwdSeq->length();) {
			MEM mem = MEM::findMEM(fwdSeq, fmidx, i, MEM::REV);

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				mem.findLocs();
				/* accept only if it is not too far */
				if(revcom_mems_pe.first.empty() || ::abs(MEM::rIndel(mtg, revcom_mems_pe.first.back(), mem)) <= maxIndelRate) {
					revcom_mems_pe.first.push_back(mem);
					i = mem.to + 1;
					revcomLoglik += mem.loglik();
				}
				else
					i++;
			}
			else
				i++;
		}
		for(int64_t i = 0; i < revSeq->length();) {
			MEM mem = MEM::findMEM(revSeq, fmidx, i, MEM::FWD);
			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = false;
			if(eval <= mem_dist(rng)) {
				mem.findLocs();
				acceptible = revcom_mems_pe.second.empty() || ::abs(MEM::rIndel(mtg, revcom_mems_pe.second.back(), mem)) <= maxIndelRate;
			}
			if(acceptible) {
					revcom_mems_pe.second.push_back(mem);
					i = mem.to + 1;
					revcomLoglik += mem.loglik();
			}
			else
				i++;
		}
	}

	if(strand & MEM::FWD != 0 && strand & MEM::REV == 0) /* fwd only */
		return sense_mems_pe;
	else if(strand & MEM::FWD == 0 && strand & MEM::REV != 0)
		return revcom_mems_pe;
	else
		return senseLoglik < revcomLoglik ? sense_mems_pe : revcom_mems_pe; /* return the most significant result */
}

} /* namespace UCSC */
} /* namespace EGriceLab */
