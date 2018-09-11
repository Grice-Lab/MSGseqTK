/*
 * MSGseqTK_main.cpp
 *  source file for MSGseqTK_main.h
 *  Created on: Sep 5, 2018
 *      Author: zhengqi
 */

#include <boost/random/uniform_01.hpp>
#include "MSGseqTK_main.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace std;

static boost::random::uniform_01<> mem_dist; /* random01 distribution for accepting MEMs */

vector<MEM> getMEMS(const PrimarySeq* seq, const FMIndex* fmidx, RNG& rng) {
	vector<MEM> mems;
	MEM mem;
	for(int64_t i = 0; i < seq->length();) {
		mem = MEM::findMEM(seq, fmidx, i);

		/* calculate MEM evalue */
		double eval = mem.evalue();
		/* accept by chance */
		bool acceptible = eval <= mem_dist(rng);
		if(acceptible) {
			mems.push_back(mem.findLocs());
			i = mem.to + 1;
		}
		else
			i++;
	}
	return mems;
}

} /* namespace EGriceLab */
} /* namespace MSGseqTK */

