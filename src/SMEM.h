/*
 * SSMEM.h
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */

#ifndef SSMEM_H_
#define SSMEM_H_

#include <vector>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <utility>
#include <algorithm>
#include "SMEM.h"
#include "PrimarySeq.h"
#include "FMDIndex.h"
#include "MetaGenome.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;

struct SMEM;
typedef vector<SMEM> SMEM_LIST;

/**
 * A Super-Maximum Exact Match class representing a
 * non-reducible match between a PrimarySeq read and a MetaGenome/FMD-Iidex target
 */
struct SMEM {
	/** default constructor */
	SMEM() = default;

	/** construct an SMEM with all info */
	SMEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size,
			const vector<GLoc>& locs)
	: seq(seq), mtg(mtg), fmdidx(fmdidx),
	  from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size), locs(locs)
	{ 	}

	/** construct an SMEM with all info but not locs */
	SMEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size)
	: seq(seq), mtg(mtg), fmdidx(fmdidx), from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size)
	{   }

	/* member methods */
	/* static member methods */
	/** test whether two MEM overlap on the seq */
	static bool isOverlap(const SMEM& lhs, const SMEM& rhs) {
		return lhs.from < rhs.to && lhs.to > rhs.from;
	}

	/** test whether two MEM is compatitable
	 * @return  true if they are ordered
	 */
	static bool isCompatitable(const SMEM& lhs, const SMEM& rhs) {
		return lhs.to <= rhs.from;
	}

	/** test whether two Loc is compatitable
	 * @return  true if they are ordered
	 */
	static bool isCompatitable(const MetaGenome* mtg, const Loc& lhs, const Loc& rhs) {
		return mtg->getLocId(lhs.start) == mtg->getLocId(rhs.start) &&
				mtg->getStrand(lhs.start) == mtg->getStrand(rhs.start) &&
				lhs.end <= rhs.start;
	}

	/** get the seq-distance of two mem,
	 * return 0 if they are overlapping
	 */
	static int64_t seqDist(const SMEM& lhs, const SMEM& rhs) {
		return isOverlap(lhs, rhs) ? 0 : lhs.from < rhs.from ? rhs.from - lhs.to + 1 : lhs.from - rhs.to + 1;
	}

	/** get the DB-distance of two mem,
	 * return 0 if they are overlapping, or SIZE_MAX no compatitable locs found
	 */
	static int64_t dbDist(const SMEM& lhs, const SMEM& rhs);

	/** get the number of indeals of two MEM
	 * return positive number if insertion, negative if deletion, or 0 if none
	 */
	static int64_t nIndel(const SMEM& lhs, const SMEM& rhs) {
		return seqDist(lhs, rhs) - dbDist(rhs, rhs);
	}

	/** get the indel rate relative to the size of their mapped seq
	 * return positive number if insertion, negative if deletion, or 0 if none
	 */
	static double rIndel(const SMEM& lhs, const SMEM& rhs) {
		return static_cast<double> (nIndel(lhs, rhs)) / lhs.seq->length();
	}

	/** get length of this SMEM */
	int64_t length() const {
		return to - from;
	}

	/** test whether this SMEM is empty */
	bool empty() const {
		return length() == 0;
	}

	/** reset this SMEM to initial state */
	void reset() {
		from = 0;
		to = 0;
		fwdStart = 0;
		revStart = 0;
		size = 0;
		logP = 0;
		seq = nullptr;
		mtg = nullptr;
		fmdidx = nullptr;
		locs.clear();
	}

	/**
	 * fill the locs information of this SMEM only on the fwd reference strand
	 * @param maxNLocs  max # of locs to find
	 */
	SMEM& findLocs(size_t maxNLocs = MAX_NLOCS);

	/**
	 * evaluate the log-probality of this SMEM
	 * @return  log-pvalue of observing this SMEM by chance, using base-frequency only
	 */
	SMEM& evaluate();

	/**
	 * get loglikelihood of this SMEM, use pre-evaluate logP value
	 * @return  log-pvalue of observing this SMEM by chance, using base-frequency only
	 */
	double loglik() const {
		return logP;
	}

	/** get the pvalue of observing this SMEM by random */
	double pvalue() const {
		return std::exp(loglik());
	}

	/** get the evalue of observing this SMEM by random */
	double evalue() const {
		return fmdidx->length() * pvalue();
	}

	/** get log-evalue of observing this SMEM by random */
	double logevalue() const {
		return std::log(fmdidx->length()) + loglik();
	}

	/** write this SMEM to text output */
	ostream& write(ostream& out) const;

	/* member fields */
	const PrimarySeq* seq = nullptr;
	const MetaGenome* mtg = nullptr;
	const FMDIndex* fmdidx = nullptr;
	int64_t from = 0; /* 0-based relative start on seq */
	int64_t to = 0;   /* 1-based relative end on seq */
	int64_t fwdStart = 0; /* 0-based start position on SA of fwd match */
	int64_t revStart = 0;   /* 0-based end position on SA of rev match */
	int64_t size = 0; /* size of match region on SA */
	double logP = 0; /* log-probability (loglik) of observing this SMEM by chance */
	vector<GLoc> locs; /* locs of mathces */

	/* static fields */
	static const size_t MAX_NLOCS = 256;
	static const double DEFAULT_MAX_EVALUE;

	/* static methods */
	/**
	 * get an SMEM_LIST of a given seq starting at given position relative to the seq by forward/backward extensions
	 * the locs will not be filled by this method
	 * @param seq  seq to search, must be in reversed orientation of this FM-index
	 * @param mtg  MetaGenome
	 * @param fmdidx  FMD-index
	 * @param smems  SMEM_LIST to store results
	 * @param from  start on the seq, will be updated after search
	 * @param to  end on seq, will be updated after search
	 * @return  SMEM_LIST found at this position
	 */
	static SMEM_LIST findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t &from, int64_t& to);

	/** evaluate an SMEM_LIST */
	static SMEM_LIST& evaluate(SMEM_LIST& smems) {
		for(SMEM& smem : smems)
			smem.evaluate();
		return smems;
	}

	/** filter an SMEM_LIST by removing SMEM with evalue great than threshold */
	static SMEM_LIST& filter(SMEM_LIST& smems, double maxEvalue = DEFAULT_MAX_EVALUE) {
		smems.erase(std::remove_if(smems.begin(), smems.end(), [=](const SMEM& smem) { return smem.evalue() > maxEvalue; }),
				smems.end());
		return smems;
	}

	/** find locs for all SMEMS */
	static SMEM_LIST& findLocs(SMEM_LIST& smems, uint32_t maxNLocs = SMEM::MAX_NLOCS) {
		for(SMEM& smem : smems)
			smem.findLocs(maxNLocs);
		return smems;
	}

	/* non-member methods */
	friend ostream& operator<<(ostream& out, const SMEM& smem);

	friend SMEM_LIST operator+(const SMEM_LIST& lhs, const SMEM_LIST& rhs);
};

inline ostream& operator<<(ostream& out, const SMEM& smem) {
	return smem.write(out); /* call virtual member method */
}

inline SMEM_LIST operator+(const SMEM_LIST& lhs, const SMEM_LIST& rhs) {
	SMEM_LIST smems_merged(lhs);
	smems_merged.insert(smems_merged.end(), rhs.begin(), rhs.end());
	return smems_merged;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SSMEM_H_ */
