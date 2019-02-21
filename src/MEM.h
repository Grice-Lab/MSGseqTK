/*
 * MEM.h
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#ifndef SRC_MEM_H_
#define SRC_MEM_H_

#include <vector>
#include <cmath>
#include <algorithm>
#include "PrimarySeq.h"
#include "FMDIndex.h"
#include "GLoc.h"
#include "MetaGenome.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {

/*
 * Maximal Exact Match betwen a PrimarySeq and an FM-index database
 */
using std::vector;

struct MEM {
	/** default constructor */
	MEM() = default;

	/** construct an MEM with all info */
	MEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size,
			const vector<GLoc>& locs)
	: seq(seq), mtg(mtg), fmdidx(fmdidx),
	  from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size), locs(locs)
	{ 	}

	/** construct an MEM with all info but not locs */
	MEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size)
	: seq(seq), mtg(mtg), fmdidx(fmdidx), from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size)
	{   }

	/* member methods */
	/** get length of this MEM */
	int64_t length() const {
		return to - from;
	}

	/** test whether this MEM is empty */
	bool empty() const {
		return length() == 0;
	}

	/** reset this MEM to initial state */
	void reset() {
		from = 0;
		to = 0;
		fwdStart = 0;
		revStart = 0;
		size = 0;
		seq = nullptr;
		mtg = nullptr;
		fmdidx = nullptr;
		locs.clear();
	}

	/**
	 * fill the locs information of this MEM only on the fwd reference strand
	 * @param maxNLocs  max # of locs to find
	 */
	MEM& findLocs(size_t maxNLocs = MAX_NLOCS);

	/**
	 * evaluate the log-probality of this MEM
	 * @return  log-pvalue of observing this MEM by chance, using base-frequency only
	 */
	MEM& evaluate();

	/**
	 * get loglikelihood of this MEM, use pre-evaluate logP value
	 * @return  log-pvalue of observing this MEM by chance, using base-frequency only
	 */
	double loglik() const {
		return logP;
	}

	/** get the pvalue of observing this MEM by random */
	double pvalue() const {
		return std::exp(loglik());
	}

	/** get the evalue of observing this MEM by random */
	double evalue() const {
		return fmdidx->length() * pvalue();
	}

	/** get log-evalue of observing this MEM by random */
	double logevalue() const {
		return std::log(fmdidx->length()) + loglik();
	}

	/** write this MEM to text output */
	ostream& write(ostream& out) const;

	/* static member methods */
	/** test whether two MEM overlap on the seq */
	static bool isOverlap(const MEM& lhs, const MEM& rhs) {
		return lhs.from < rhs.to && lhs.to > rhs.from;
	}

	/** test whether two MEM is compatitable
	 * @return  true if they are ordered
	 */
	static bool isCompatitable(const MEM& lhs, const MEM& rhs) {
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
	static int64_t seqDist(const MEM& lhs, const MEM& rhs) {
		return isOverlap(lhs, rhs) ? 0 : lhs.from < rhs.from ? rhs.from - lhs.to + 1 : lhs.from - rhs.to + 1;
	}

	/** get the DB-distance of two mem,
	 * return 0 if they are overlapping, or SIZE_MAX no compatitable locs found
	 */
	static int64_t dbDist(const MEM& lhs, const MEM& rhs);

	/** get the number of indeals of two MEM
	 * return positive number if insertion, negative if deletion, or 0 if none
	 */
	static int64_t nIndel(const MEM& lhs, const MEM& rhs) {
		return seqDist(lhs, rhs) - dbDist(rhs, rhs);
	}

	/** get the indel rate relative to the size of their mapped seq
	 * return positive number if insertion, negative if deletion, or 0 if none
	 */
	static double rIndel(const MEM& lhs, const MEM& rhs) {
		return static_cast<double> (nIndel(lhs, rhs)) / lhs.seq->length();
	}

	/* member fields */
	const PrimarySeq* seq = nullptr;
	const MetaGenome* mtg = nullptr;
	const FMDIndex* fmdidx = nullptr;
	int64_t from = 0; /* 0-based relative start on seq */
	int64_t to = 0;   /* 1-based relative end on seq */
	int64_t fwdStart = 0; /* 0-based start position on SA of fwd match */
	int64_t revStart = 0;   /* 0-based end position on SA of rev match */
	int64_t size = 0; /* size of match region on SA */
	double logP = 0; /* log-probability (loglik) of observing this MEM by chance */
	vector<GLoc> locs; /* locs of mathces */

	/* static fields */
	static const size_t MAX_NLOCS = 256;

	/* static methods */
	/**
	 * get an MEM of a given seq starting at given position relative to the seq
	 * the locs will not be filled by this method
	 * @param seq  seq to search, must be in reversed orientation of this FM-index
	 * @param i  relative position of the seq
	 * @param strand  which direction/strand to search
	 */
	static MEM findMEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from = 0, int64_t to = INT64_MAX, GLoc::STRAND dir = GLoc::FWD) {
		return dir == GLoc::FWD ? findMEMfwd(seq, mtg, fmdidx, from) :
				findMEMrev(seq, mtg, fmdidx, to);
	}

	static MEM findMEMfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx, int64_t from = 0);

	static MEM findMEMrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx, int64_t to = INT64_MAX);

	/* non-member methods */
	friend ostream& operator<<(ostream& out, const MEM& mem);
};

inline ostream& operator<<(ostream& out, const MEM& mem) {
	return mem.write(out); /* call virtual member method */
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_MEM_H_ */
