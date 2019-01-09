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
#include <algorithm>
#include "MEM.h"
#include "PrimarySeq.h"
#include "FMIndex.h"
#include "MetaGenome.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::pair;
typedef boost::random::mt11213b RNG;

struct MEMS;
typedef pair<MEMS, MEMS> MEMS_PE;

struct MEMS : public vector<MEM> {
public:
	/* constructors */
	/** default constructor */
	MEMS() = default;

	/* member methods */
	/** evaluate each MEM in this MEMS */
	MEMS& evaluate();

	/** find locs for all MEM in this MEMS */
	MEMS& findLocs(size_t maxNLocs = MEM::MAX_NLOCS);

	/** get the loglik of the entire MEMS if it is from random matches */
	double loglik() const;

	/** get the probability of MEMS if it is a correct (non-random) match */
	double Pr() const {
		return 1 - ::exp(loglik());
	}

	/** get the total length of this MEMS */
	uint64_t length() const;

	/**
	 * locate the best MEM index with lowest pvalue or likelihood
	 * @return  the MEM with smallest loglik, or -1 if empty
	 */
	size_t bestMEMIndex() const;

	/** locate the best MEM in this MEMS */
	const MEM& bestMEM() const {
		return (*this)[bestMEMIndex()];
	}

	/** write this MEMS to text output */
	ostream& write(ostream& out) const;

	/** get strand of this MEMS */
	MEM::STRAND getStrand() const {
		return front().strand;
	}

	/* static methods */
	/**
	 * get one MEMS by MCMC sampling matches between db and seq
	 * @param seq  primary sequence to search
	 * @param fmidx  FM-index
	 * @param rng  random-number generator
	 * @param strand  searching strand
	 * @return  a vector of ordered MEMs
	 */
	static MEMS sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
			RNG& rng, MEM::STRAND strand);

	/**
	 * get best MEMS by trying different strands and multiple random seeds
	 * @param strand  different strands, 1 for FWD, 2 for REV, 3 for both
	 * @param nSeed  # of different seeds
	 */
	static MEMS getBestMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
			RNG& rng, int strand, int nSeed = DEFAULT_NSEED);

	/**
	 * get best MEMS_PE by trying different strands and multiple random seeds
	 * @param strand  different strands, 1 for FWD, 2 for REV, 3 for both
	 * @param nSeed  # of different seeds
	 */
	static MEMS_PE getBestMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx,
			RNG& rng, int strand, int nSeed = DEFAULT_NSEED);

	/** find locs for MEMS_PE */
	static MEMS_PE& findLocs(MEMS_PE& mems_pe, size_t maxNLocs = MEM::MAX_NLOCS);

//	/** filter MEMs by removing MEM with exceeding number of locs */
//	MEMS& filterLocs(size_t maxNLocs = MAX_MEM_NLOCS) {
//		erase(std::remove_if(begin(), end(), [=] (const MEM& mem) { return mem.locs.size() > maxNLocs; }));
//		return *this;
//	}

	/**
	 * filter MEMS by removing incompatitable MEM's locs that is not on the same genome, or same chromosome and with not too much indels
	 * this algorithm starts from the best MEM (lowest loglik), then sweep to both ends
	 * @param mtg  MetaGenome data
	 * @param maxIndel  max indel-rate allowed
	 */
	MEMS& filterLocs(const MetaGenome& mtg, double maxIndel);

	/**
	 * get the loglik of a paired-end MEMS
	 */
	static double loglik(const MEMS_PE& mems_pe) {
		return mems_pe.first.loglik() + mems_pe.second.loglik();
	}

	/** get strand of a MEMS_PE */
	static MEM::STRAND getStrand(const MEMS_PE& mems_pe) {
		return mems_pe.first.getStrand();
	}

	/* static fields */
	static boost::random::uniform_01<> mem_dist; /* random01 distribution for accepting MEMs */
	static const size_t MAX_MEM_NLOCS = 4;
	static const int DEFAULT_NSEED = 1;

	/* non-member methods */
	friend ostream& operator<<(ostream& out, const MEMS& mems);
};

inline MEMS& MEMS::evaluate() {
	for(MEM& mem : *this)
		mem.evaluate();
	return *this;
}

inline double MEMS::loglik() const {
	double logP = 0;
	for(const MEM& mem : *this)
		logP += mem.loglik();
	return logP;
}

inline ostream& operator<<(ostream& out, const MEMS& mems) {
	return mems.write(out); /* call virtual member method */
}

inline MEMS& MEMS::findLocs(size_t maxNLocs) {
	for(MEM& mem : *this)
		mem.findLocs(maxNLocs);
	return *this;
}

inline MEMS_PE& MEMS::findLocs(MEMS_PE& mems_pe, size_t maxNLocs) {
	mems_pe.first.findLocs(maxNLocs);
	mems_pe.second.findLocs(maxNLocs);
	return mems_pe;
}

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* SRC_MEMS_H_ */
