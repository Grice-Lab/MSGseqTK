/*
 * MEM.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#include <cmath>
#include "MEM.h"
#include "Stats.h"

namespace EGriceLab {
namespace MSGseqClean {

using namespace EGriceLab::Math;

double MEM::loglik(const Vector4d& baseFreq) const {
	Vector4d logBaseFreq = baseFreq.array().log();
	double ll = 0; /* ll is always in base 10 in phred system */
	for(uint64_t i = from; i < to; ++i) {
		ll += logBaseFreq((*seq)[i] - DNAalphabet::A); /* use observed base frequency */

		if(qual != nullptr) /* quality is in use */
			ll += ::log(1 - q2p((*qual)[i])); /* matching liklihood is 1 - error */
	}
	return ll;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */


