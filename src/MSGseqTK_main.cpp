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

vector<MEM> getMEMS(const FMIndex& fmidx, const PrimarySeq& seq, RNG& rng, bool ignoreQual) {
	vector<MEM> mems;
	MEM mem;
	for(int64_t i = 0; i < seq.length();) {
		if(!ignoreQual && seq.hasQual())
			mem = fmidx.findMEM(seq.getSeq(), seq.getQual(), i);
		else
			mem = fmidx.findMEM(seq.getSeq(), i);

		/* calculate MEM evalue */
		double eval = mem.evalue();
		/* accept by chance */
		bool acceptible = mem_dist(rng) <= eval;
		if(acceptible) {
			mems.push_back(mem);
			i = mem.to + 1;
		}
		else
			i++;
	}
	return mems;
}

} /* namespace EGriceLab */
} /* namespace MSGseqTK */

