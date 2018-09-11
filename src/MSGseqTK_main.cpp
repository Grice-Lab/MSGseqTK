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
	if(strand | MEM::FWD != 0) {
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
	if(strand | MEM::REV != 0) {
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

	if(strand | MEM::FWD != 0 && strand | MEM::REV == 0) /* fwd only */
		return fwdMems;
	else if(strand | MEM::FWD == 0 && strand | MEM::REV != 0)
		return revMems;
	else
		return fwdLoglik > revLoglik ? fwdMems : revMems;
}

} /* namespace EGriceLab */
} /* namespace MSGseqTK */

