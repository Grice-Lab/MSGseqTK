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
#include "SeedPair.h"
#include "PrimarySeq.h"
#include "FMDIndex.h"
#include "MetaGenome.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::pair;

/**
 * A Super-Maximum Exact Match class representing a
 * non-reducible match between a PrimarySeq read and a MetaGenome/FMD-Iidex target
 */
class SMEM {
public:
	/** default constructor */
	SMEM() = default;

	/** construct an SMEM with all info */
	SMEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size)
	: seq(seq), mtg(mtg), fmdidx(fmdidx),
	  from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size)
	{
		evaluate();
	}

	/* member methods */
	/** getters and setters */
	const PrimarySeq* getSeq() const {
		return seq;
	}

	const FMDIndex* getFmdidx() const {
		return fmdidx;
	}

	const MetaGenome* getMtg() const {
		return mtg;
	}

	int64_t getFrom() const {
		return from;
	}

	int64_t getTo() const {
		return to;
	}

	int64_t getSize() const {
		return size;
	}

	double loglik() const {
		return logP;
	}

	double pvalue() const {
		return std::exp(loglik());
	}

	double evalue() const {
		return fmdidx != nullptr ? fmdidx->length() * pvalue() : 0;
	}

	/** get length of this SMEM */
	int64_t length() const {
		return to - from;
	}

	/** test whether this SMEM is empty */
	bool empty() const {
		return length() == 0;
	}

	/**
	 * get SeedPairs of this SMEM
	 * @return  mapped seeds that always on FWD strand of target
	 */
	SeedList getSeeds() const;

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

	/* static methods */
	/**
	 * find one longest SMEM of a given seq starting at given position relative to seq by forward or backward extension
	 */
	static SMEM findSMEM(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t& from, int64_t& to, GLoc::STRAND dir = GLoc::FWD) {
		return dir == GLoc::FWD ? findSMEMfwd(seq, mtg, fmdidx, from, to) :
				findSMEMrev(seq, mtg, fmdidx, from, to);
	}

	/**
	 * find one longest SMEM of a given seq starting at given position relative to seq by forward extension
	 */
	static SMEM findSMEMfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t& to);

	/**
	 * find one longest SMEM of a given seq starting at given position relative to seq by backward extension
	 */
	static SMEM findSMEMrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t& from, int64_t to);

private:
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

public:
	/* static fields */
	static const size_t MAX_NSEEDS = 256;
	static const size_t MAX_NCHAINS = 256;

	friend class SMEMS;
};

class SMEMS;
typedef pair<SMEMS, SMEMS> SMEMS_PE;

/**
 * a SMEMS is a list SMEM with additional methods
 */
class SMEMS : public vector<SMEM> {
public:
	/* member methods */
	/** get total loglik */
	double loglik() const;

	/** get total pvalue */
	double pvalue() const {
		return std::exp(loglik());
	}

	/** get total evalue */
	double evalue() const {
		return front().fmdidx->length() * pvalue();
	}

	/** get best (min) loglik */
	double bestLoglik() const;

	/** get best (min) pvalue */
	double bestPvalue() const {
		return std::exp(bestLoglik());
	}

	/** get best (min) evalue */
	double bestEvalue() const {
		return front().fmdidx->length() * bestPvalue();
	}

	/**
	 * filter this SMEMS list by evalue
	 */
	SMEMS& filter(double maxEvalue = DEFAULT_MAX_EVALUE) {
		erase(std::remove_if(begin(), end(),
				[=](const SMEM& mem) { return mem.evalue() > maxEvalue; }),
				end());
		return *this;
	}

	/* static methods */
	/**
	 * find longest SMEMS of a given seq using step-wise forward/backward searches with required evalue threshold
	 */
	static SMEMS findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE, GLoc::STRAND dir = GLoc::FWD) {
		return dir == GLoc::FWD ? findSMEMSfwd(seq, mtg, fmdidx, maxEvalue) :
				findSMEMSrev(seq, mtg, fmdidx, maxEvalue);
	}

	/**
	 * find longest SMEMS of a given seq using step-wise backward searches with required evalue threshold
	 */
	static SMEMS findSMEMSfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	/**
	 * find longest SMEMS of a given seq using step-wise forward searches with required evalue threshold
	 */
	static SMEMS findSMEMSrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	/**
	 * find all SMEMS of a given seq starting at given position relative to the seq by forward/backward extensions
	 */
	static SMEMS findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t& from, int64_t& to);

	/**
	 * find all SMEMS of a given seq using step-wise forward/backward searches
	 */
	static SMEMS findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	/**
	 * get a SeedList of a given seq using step-wise forward/backward searches
	 * seeds will be filtered and sorted
	 */
	static SeedList findSeeds(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			double maxEvalue = DEFAULT_MAX_EVALUE);

	/**
	 * find SMEMS_PE for paired-end reads
	 */
	static SMEMS_PE findSMEMS_PE(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq,
			const MetaGenome* mtg, const FMDIndex* fmdidx, double maxEvalue = DEFAULT_MAX_EVALUE) {
		return SMEMS_PE(findSMEMS(fwdSeq, mtg, fmdidx, maxEvalue), findSMEMS(revSeq, mtg, fmdidx, maxEvalue));
	}

	/**
	 * find SeedListPE for pair-end reads
	 */
	static SeedListPE findSeedsPE(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq,
			 const MetaGenome* mtg, const FMDIndex* fmdidx, double maxEvalue = DEFAULT_MAX_EVALUE) {
		return SeedListPE(findSeeds(fwdSeq, mtg, fmdidx, maxEvalue), findSeeds(revSeq, mtg, fmdidx, maxEvalue));
	}

	/** get loglik for SMEMS_PE */
	static double loglik(const SMEMS_PE& smemsPE) {
		return smemsPE.first.loglik() + smemsPE.second.loglik();
	}

	/* static fields */
	static const double DEFAULT_MAX_EVALUE;
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SSMEM_H_ */
