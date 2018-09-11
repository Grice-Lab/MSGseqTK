/*
 * MSGseqTK_main.h
 *  Declaration of MSGseqTK core algorithms
 *
 *  Created on: Sep 5, 2018
 *      Author: zhengqi
 *      Since: v1.1
 */

#ifndef SRC_MSGSEQTK_MAIN_H_
#define SRC_MSGSEQTK_MAIN_H_

#include <string>
#include <vector>
#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include "MSGseqTK.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::vector;

typedef boost::random::mt11213b RNG;

/**
 * get MEMS by MCMC sampling matches between db and seq
 * @param fmidx  FM-index
 * @param seq  primary sequence to search
 * @param ignoreQual  whether to force ignoring quality info even if it is presented in the sequence
 */
vector<MEM> getMEMS(const PrimarySeq* seq, const FMIndex* fmidx, RNG& rng);

/**
 * filter MEMs by removing incompatitable MEMs that is not on the same genome, same chromosome and with not too much indels
 * @param mems  vector of MEMs ordered by their relative location on seq
 * @param mtg  MetaGenome data
 * @param maxIndel  max indel-rate allowed
 */
vector<MEM>& filterMEMS(vector<MEM>& mems, const MetaGenome& mtg, double maxIndel);


} /* namespace EGriceLab */
} /* namespace MSGseqTK */

#endif /* SRC_MSGSEQTK_MAIN_H_ */
