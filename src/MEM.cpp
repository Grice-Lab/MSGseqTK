/*
 * MEM.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#include <cmath>
#include <cassert>
#include "MEM.h"
#include "Stats.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace EGriceLab::Math;

double MEM::logP() const {
	double ll = 0; /* ll is always in base 10 in phred system */
	for(uint64_t i = from; i < to; ++i) {
		ll += ::log(B[(*seq)[i]]) - ::log(N); /* use observed base frequency */

		if(qual != nullptr) /* quality is in use */
			ll += ::log(1 - q2p((*qual)[i])); /* matching liklihood is 1 - error */
	}
	return ll;
}

int64_t MEM::readDist(const MEM& mem1, const MEM& mem2) {
	assert(mem1.seq == mem2.seq);
	if(isOverlap(mem1, mem2))
		return 0;
	else
		return mem1.from < mem2.from ? mem2.from - mem1.to + 1 : mem1.from - mem2.to + 1;
}

int64_t MEM::dbDist(const MEM& mem1, const MEM& mem2) {
	int64_t minD = -1;
	for(const Loc& loc1 : mem1.locs) {
		for(const Loc& loc2 : mem2.locs) {
			int64_t d = Loc::dist(loc1, loc2);
			if(minD == -1 || d < minD)
				minD = d;
		}
	}
	return minD;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */


