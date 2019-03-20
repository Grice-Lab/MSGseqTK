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
#include <cassert>
#include <utility>
#include <algorithm>
#include "MEM.h"
#include "PrimarySeq.h"
#include "FMDIndex.h"
#include "MetaGenome.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;

struct MEMS : public vector<MEM> {
	/* constructors */
	/** default constructor */
	MEMS() = default;

	/** destructor */
	virtual ~MEMS() {  }

	/* member methods */
	/** get seq */
	const PrimarySeq* getSeq() const {
		return front().seq;
	}

	/** get FMD-index */
	const FMDIndex* getFMDIndex() const {
		return front().fmdidx;
	}

	/** evaluate each MEM in this MEMS */
	MEMS& evaluate();

	/** find locs for all MEMs */
	MEMS& findLocs(size_t maxNLocs = MEM::MAX_NLOCS);

	/** get the loglik of the entire MEMS if it is from random matches */
	double loglik() const;

	/** get the evalue of observing this MEMS by random */
	double evalue() const {
		return getFMDIndex()->length() * std::exp(loglik());
	}

	/** get the log-evalue of observing this MEMS by random */
	double logevalue() const {
		return std::log(getFMDIndex()->length()) + loglik();
	}

	/** get the total length of this MEMS */
	int64_t length() const;

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

	/** append another MEMS to this MEMS */
	MEMS& operator+=(const MEMS& other);

	/* static methods */
	/**
	 * search MEMS for given read
	 * @param seq  read to search
	 * @param fmdidx  FMD-index
	 * @param maxEvalue  maxEvalue criteria for a significant MEM
	 * @param dir  searching direction
	 * @return  an ordered MEMs
	 */
	static MEMS searchMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE, GLoc::STRAND dir = GLoc::FWD) {
		return dir == GLoc::FWD ? searchMEMSfwd(seq, mtg, fmdidx, maxEvalue) :
				searchMEMSrev(seq, mtg, fmdidx, maxEvalue);
	}

	static MEMS searchMEMSfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	static MEMS searchMEMSrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	/* static fields */
	static const size_t MAX_MEM_NLOCS = 256;
	static const double DEFAULT_MAX_EVALUE;

	/* non-member methods */
	/** formatted output */
	friend ostream& operator<<(ostream& out, const MEMS& mems);
	/** merge two MEMS */
	friend MEMS operator+(const MEMS& lhs, const MEMS& rhs);
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

	/* static methods */
	/**
	 * sample MEMS_PE for pair-end reads
	 * @param fwdSeq  forward read
	 * @param revSeq  reverse read
	 * @param fmdidx  FMD-index
	 * @param maxEvalue  maxEvalue criteria for a significant MEM
	 * @param dir  direction for sampling
	 */
	static MEMS_PE sampleMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq,
			 const MetaGenome* mtg, const FMDIndex* fmdidx,
			 double maxEvalue = MEMS::DEFAULT_MAX_EVALUE, GLoc::STRAND dir = GLoc::FWD) {
		return dir == GLoc::FWD ?
				MEMS_PE(MEMS::searchMEMSfwd(fwdSeq, mtg, fmdidx, maxEvalue), MEMS::searchMEMSfwd(revSeq, mtg, fmdidx, maxEvalue)) :
				MEMS_PE(MEMS::searchMEMSrev(fwdSeq, mtg, fmdidx, maxEvalue), MEMS::searchMEMSrev(revSeq, mtg, fmdidx, maxEvalue));
	}

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

inline MEMS operator+(const MEMS& lhs, const MEMS& rhs) {
	MEMS mems(lhs);
	mems += rhs;
	return mems;
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
