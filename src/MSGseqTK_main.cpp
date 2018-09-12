/*
 * MSGseqTK_main.cpp
 *  source file for MSGseqTK_main.h
 *  Created on: Sep 5, 2018
 *      Author: zhengqi
 */

#include <boost/random/uniform_01.hpp>
#include <cassert>
#include "MSGseqTK_main.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace std;

static boost::random::uniform_01<> mem_dist; /* random01 distribution for accepting MEMs */

vector<MEM> getMEMS(const PrimarySeq* seq, const FMIndex* fmidx, RNG& rng, int strand) {
	assert(strand != 0);
	vector<MEM> fwdMems, revMems;
	double fwdLoglik = 0;
	double revLoglik = 0;
	/* scan fwd strand */
	if(strand & MEM::FWD != 0) {
		for(int64_t i = 0; i < seq->length();) {
			MEM mem = MEM::findMEM(seq, fmidx, i, MEM::FWD);

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				fwdMems.push_back(mem.findLocs());
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
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				revMems.push_back(mem.findLocs());
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

MEMS_PE getMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx, RNG& rng, int strand) {
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
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				sense_mems_pe.first.push_back(mem.findLocs());
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
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				sense_mems_pe.second.push_back(mem.findLocs());
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
				revcom_mems_pe.first.push_back(mem.findLocs());
				i = mem.to + 1;
				revcomLoglik += mem.loglik();
			}
			else
				i++;
		}
		for(int64_t i = 0; i < revSeq->length();) {
			MEM mem = MEM::findMEM(revSeq, fmidx, i, MEM::FWD);

			/* calculate MEM evalue */
			double eval = mem.evalue();
			/* accept by chance */
			bool acceptible = eval <= mem_dist(rng);
			if(acceptible) {
				revcom_mems_pe.second.push_back(mem.findLocs());
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

} /* namespace EGriceLab */
} /* namespace MSGseqTK */

