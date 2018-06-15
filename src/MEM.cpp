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

double MEM::loglik() const {
	double ll = 0; /* ll is always in base 10 in phred system */
	for(uint64_t i = from; i < to; ++i) {
		uint8_t q = qual != nullptr ? (*qual)[i] : QualStr::DEFAULT_Q_SCORE;
		ll += ::log(1 - q2p(q)); /* matching liklihood is 1 - error */
	}
	return ll;
}

Vector4d MEM::baseFreq() const {
	Vector4d baseFreq = Vector4d::Zero();
	for(DNAseq::value_type b : *seq)
		if(DNAalphabet::isBase(b))
			baseFreq(b - DNAalphabet::A)++;
	return baseFreq / baseFreq.sum();
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */


