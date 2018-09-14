/*
 * MEMS.h
 *  Minimal Exact Match Seeds
 *  Created on: Sep 12, 2018
 *      Author: zhengqi
 */

#ifndef SRC_MEMS_H_
#define SRC_MEMS_H_

#include <vector>
#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <cassert>
#include <utility>
#include "MEM.h"
#include "PrimarySeq.h"
#include "FMIndex.h"
#include "MetaGenome.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::pair;
typedef boost::random::mt11213b RNG;

class MEMS;
typedef pair<MEMS, MEMS> MEMS_PE;

class MEMS : public vector<MEM> {
public:
	/* constructors */
	/** default constructor */
	MEMS() = default;

	/* member methods */
	/** get the loglik of the entire MEMS if it is from random matches */
	double loglik() const;

	/** get the probability of MEMS if it is a correct (non-random) match */
	double Pr() const {
		return 1 - ::exp(loglik());
	}

	/**
	 * locate the best MEM index with lowest pvalue or likelihood
	 * @return  the MEM with smallest loglik, or -1 if empty
	 */
	size_t bestMEMIndex() const;

	/** locate the best MEM in this MEMS */
	const MEM& bestMEM() const {
		return (*this)[bestMEMIndex()];
	}

	/**
	 * filter MEMS by removing incompatitable MEM's locs that is not on the same genome, or same chromosome and with not too much indels
	 * this algorithm starts from the best MEM (lowest loglik), then sweep to both ends
	 * @param mtg  MetaGenome data
	 * @param maxIndel  max indel-rate allowed
	 */
	MEMS& filterLocs(const MetaGenome& mtg, double maxIndel);

	/* static methods */
	/**
	 * get MEMS by MCMC sampling matches between db and seq
	 * @param seq  primary sequence to search
	 * @param fmidx  FM-index
	 * @param mtg  MetaGenome
	 * @param rng  random-number generator
	 * @return  a vector of MEMs ordered by from
	 */
	static MEMS sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx, const MetaGenome* mtg,
			RNG& rng, int strand, double maxIndelRate);

	/**
	 * get paired-end MEMS by MCMC sampling matches between db and seq
	 * @param fwdSeq  fwd-seq to search
	 * @param revSeq  rev-seq to search
	 * @param fmidx  FM-index
	 * @param rng  random-number generator
	 * @return  an MEM_PE of MEMs ordered by from
	 */
	MEMS_PE sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx, const MetaGenome* mtg,
			RNG& rng, int strand, double maxIndelRate);

	/* static fields */
	static boost::random::uniform_01<> mem_dist; /* random01 distribution for accepting MEMs */
};

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* SRC_MEMS_H_ */
