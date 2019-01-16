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

struct MEMS : public vector<MEM> {
public:
	/* constructors */
	/** default constructor */
	MEMS() = default;

	/* member methods */
	/** get seq */
	const PrimarySeq* getSeq() const {
		return front().seq;
	}

	/** get FM-index */
	const FMIndex* getFMIndex() const {
		return front().fmidx;
	}

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
	 */
	static MEMS sampleMEMS(const PrimarySeq* seq, const FMIndex* fmidx,
			RNG& rng, int strand);


	/* static fields */
	static boost::random::uniform_01<> mem_dist; /* random01 distribution for accepting MEMs */
	static const size_t MAX_MEM_NLOCS = 4;

	/* non-member methods */
	friend ostream& operator<<(ostream& out, const MEMS& mems);
};

/** convenient wrapper class for MEMS paired-end (PE) */
struct MEMS_PE {
	/* constructors */
	/** default constructor */
	MEMS_PE() = default;

	/** construct MEMS_PE given both pair MEMS */
	MEMS_PE(const MEMS& fwdMems, const MEMS& revMems) : fwdMems(fwdMems), revMems(revMems)
	{  }

	/* member methods */
	/** find locs for both fwd and rev MEMS */
	MEMS_PE& findLocs(size_t maxNLocs = MEM::MAX_NLOCS);

	/** get the loglik of this MEMS_PE */
	double loglik() const {
		return fwdMems.loglik() + revMems.loglik();
	}

	/** get forward strand of this MEMS_PE */
	MEM::STRAND getFwdStrand() const {
		return fwdMems.getStrand();
	}

	/** get reverse strand of this MEMS_PE */
	MEM::STRAND getRevStrand() const {
		return revMems.getStrand();
	}

	/** get strand of this MEMS_PE, alias to getFwdStrand */
	MEM::STRAND getStrand() const {
		return getFwdStrand();
	}

	/* static methods */
	/**
	 * get best MEMS_PE by trying different strands and multiple random seeds
	 * @param strand  different strands, 1 for FWD, 2 for REV, 3 for both
	 */
	static MEMS_PE sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq, const FMIndex* fmidx,
			RNG& rng, int strand);

	/* member fields */
	MEMS fwdMems;
	MEMS revMems;
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

inline MEMS_PE& MEMS_PE::findLocs(size_t maxNLocs) {
	fwdMems.findLocs(maxNLocs);
	revMems.findLocs(maxNLocs);
	return *this;
}

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* SRC_MEMS_H_ */
