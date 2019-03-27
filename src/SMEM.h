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
using std::pair;

struct SMEM;
typedef vector<SMEM> SMEMS;
typedef pair<SMEMS, SMEMS> SMEMS_PE;

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
	{
		evaluate();
	}

	/** construct an SMEM with all info but not locs */
	SMEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size)
	: seq(seq), mtg(mtg), fmdidx(fmdidx), from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size)
	{
		evaluate();
	}

	/* member methods */
	/* static member methods */
	/** test whether two MEM overlap on the seq */
	static bool isOverlap(const SMEM& lhs, const SMEM& rhs) {
		return lhs.from < rhs.to && lhs.to > rhs.from;
	}

	/** test whether two MEM is compatitable
	 * @return  true if they are linearly compatitable
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

	/** test whether one SMEM contains the other */
	static bool contains(const SMEM& lhs, const SMEM& rhs) {
		return lhs.from <= rhs.from && lhs.to >= rhs.to;
	}

	/** testh wheter one SMEM is contained in the other */
	static bool contained(const SMEM& lhs, const SMEM& rhs) {
		return contains(rhs, lhs);
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
	 * forward extend this SMEM at current end
	 * @return  updated SMEM
	 */
	SMEM& fwdExt() {
		fmdidx->fwdExt(fwdStart, revStart, size, seq->getBase(to));
		fwdEvaluate();
		to++;
		return *this;
	}

	/**
	 * get a copy of forward extension of this SMEM
	 * @return  new SMEM with extended values
	 */
	SMEM fwdExt() const {
		SMEM smemExt(*this);
		return smemExt.fwdExt();
	}

	/**
	 * backward extend this SMEM at current start
	 * @return  updated SMEM
	 */
	SMEM& backExt() {
		fmdidx->backExt(fwdStart, revStart, size, seq->getBase(from - 1));
		backEvaluate();
		from--;
		return *this;
	}

	/**
	 * get a copy of backward extension of this SMEM
	 * @return  new SMEM with extended values
	 */
	SMEM backExt() const {
		SMEM smemExt(*this);
		return smemExt.backExt();
	}

	/**
	 * evaluate the log-probality of this SMEM
	 * @return  log-pvalue of observing this SMEM by chance, using base-frequency only
	 */
	SMEM& evaluate();

protected:
	/** forward evaluation during forward extension */
	SMEM& fwdEvaluate() {
		logP += fmdidx->loglik(seq->getBase(to));
		return *this;
	}

	/** backward evaluation during backward extension */
	SMEM& backEvaluate() {
		logP += fmdidx->loglik(seq->getBase(from - 1));
		return *this;
	}

public:
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
	static const size_t MAX_NCHAINS = 256;
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
	static SMEMS findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t &from, int64_t& to);

	/**
	 * get an aggregated SMEM_LIST of a given seq using step-wise forward/backward searches
	 * returned seq will be sorted based on the SMEM's lexical order (of from and to)
	 */
	static SMEMS searchSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	/** evaluate an SMEM_LIST */
	static SMEMS& evaluate(SMEMS& smems) {
		for(SMEM& smem : smems)
			smem.evaluate();
		return smems;
	}

	/** filter an SMEM_LIST by removing SMEM with evalue great than threshold */
	static SMEMS& filter(SMEMS& smems, double maxEvalue = DEFAULT_MAX_EVALUE) {
		smems.erase(std::remove_if(smems.begin(), smems.end(), [=](const SMEM& smem) { return smem.evalue() > maxEvalue; }),
				smems.end());
		return smems;
	}

	/** find locs for all SMEMS */
	static SMEMS& findLocs(SMEMS& smems, uint32_t maxNLocs = SMEM::MAX_NLOCS) {
		for(SMEM& smem : smems)
			smem.findLocs(maxNLocs);
		return smems;
	}

	/** order SMEM in lexicographical order */
	static SMEMS& sort(SMEMS& smems) {
		std::sort(smems.begin(), smems.end());
		return smems;
	}

	/** get the best (min) loglik of an SEMEM LIST */
	static double bestLoglik(const SMEMS& smems);

	/** get the best (min) evalue of an SMEM_LIST */
	static double bestEvalue(const SMEMS& smems);

	/** test whether an SMEM_LIST is compatitable */
	static bool isCompatitable(const SMEMS& smems);

	/** get the seq of an SMEM list */
	static const PrimarySeq* getSeq(const SMEMS& smems) {
		return smems.front().seq;
	}

	/** get the FMD-index of an SMEM list */
	static const FMDIndex* getFMDIndex(const SMEMS& smems) {
		return smems.front().fmdidx;
	}

	/** get the MetaGenome of an SMEM list */
	static const MetaGenome* getMetaGenome(const SMEMS& smems) {
		return smems.front().mtg;
	}

	/** get from of an ordered SMEM list */
	static int64_t getFrom(const SMEMS& smems) {
		return smems.front().from;
	}

	/** get the loglik of an SMEM list */
	static double loglik(const SMEMS& smems);

	/** get the evalue of an SMEM list */
	static double evalue(const SMEMS& smems) {
		return getFMDIndex(smems)->length() * std::exp(loglik(smems));
	}

	/** get the log-evalue of an MEMS list */
	static double logevalue(const SMEMS& smems) {
		return std::log(getFMDIndex(smems)->length()) + loglik(smems);
	}

	/** get the total length of an MEMS list */
	static int64_t length(const SMEMS& smems);

	/** write an SMEM list to text output */
	static ostream& write(const SMEMS& smems, ostream& out);

	/** find locs for an SMEM_PE */
	static SMEMS_PE& findLocs(SMEMS_PE& smemsPE, size_t maxNLocs = SMEM::MAX_NLOCS) {
		SMEM::findLocs(smemsPE.first, maxNLocs);
		SMEM::findLocs(smemsPE.second, maxNLocs);
		return smemsPE;
	}

	/** get the loglik of an SMEMS_PE */
	static double loglik(const SMEMS_PE& smemsPE) {
		return SMEM::loglik(smemsPE.first) + SMEM::loglik(smemsPE.second);
	}

	/**
	 * search SMEMS_PE for pair-end reads
	 */
	static SMEMS_PE searchSMEMS(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq,
			 const MetaGenome* mtg, const FMDIndex* fmdidx, double maxEvalue = SMEM::DEFAULT_MAX_EVALUE) {
		return SMEMS_PE(SMEM::searchSMEMS(fwdSeq, mtg, fmdidx, maxEvalue), SMEM::searchSMEMS(revSeq, mtg, fmdidx, maxEvalue));
	}

	/* non-member methods */
	/** formatted output for SMEM */
	friend ostream& operator<<(ostream& out, const SMEM& smem);
	/** formatted output for SMEM chain */
	friend ostream& operator<<(ostream& out, const SMEMS& smems);
	/** merge two SMEM lists */
	friend SMEMS operator+(const SMEMS& lhs, const SMEMS& rhs);

	/* relationship operators, all comparions are only based on from and to */
	friend bool operator<(const SMEM& lhs, const SMEM& rhs);
	friend bool operator==(const SMEM& lhs, const SMEM& rhs);

};

inline ostream& operator<<(ostream& out, const SMEM& smem) {
	return smem.write(out); /* call virtual member method */
}

inline ostream& operator<<(ostream& out, const SMEMS& smems) {
	return SMEM::write(smems, out);
}

inline SMEMS operator+(const SMEMS& lhs, const SMEMS& rhs) {
	SMEMS smems_merged(lhs);
	smems_merged.insert(smems_merged.end(), rhs.begin(), rhs.end());
	return smems_merged;
}

inline bool operator<(const SMEM& lhs, const SMEM& rhs) {
	assert(lhs.seq == rhs.seq && lhs.fmdidx == rhs.fmdidx && lhs.mtg == rhs.mtg);
	return lhs.from != rhs.from ? lhs.from < rhs.from : lhs.to < rhs.to;
}

inline bool operator==(const SMEM& lhs, const SMEM& rhs) {
	assert(lhs.seq == rhs.seq && lhs.fmdidx == rhs.fmdidx && lhs.mtg == rhs.mtg);
	return lhs.from == rhs.from && lhs.to == rhs.to;
}

inline bool operator<=(const SMEM& lhs, const SMEM& rhs) {
	return lhs < rhs || lhs == rhs;
}

inline bool operator>(const SMEM& lhs, const SMEM& rhs) {
	return rhs < lhs;
}

inline bool operator>=(const SMEM& lhs, const SMEM& rhs) {
	return !(lhs < rhs);
}

inline bool operator!=(const SMEM& lhs, const SMEM& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

/** customized hash function in std namespace */
namespace std {
template<>
class hash<EGriceLab::MSGseqTK::SMEM> {
public:
	size_t operator() (const EGriceLab::MSGseqTK::SMEM& smem) const {
		return hash(smem.from) ^ hash(smem.to);
	}
};
} /* namespace std */

#endif /* SSMEM_H_ */
