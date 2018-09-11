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
#include "Loc.h"
#include "PrimarySeq.h"
#include "FMIndex.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {

/*
 * Maximal Exact Match betwen a PrimarySeq and an FM-index database
 */
using std::vector;

struct MEM {
	/* nested types and enums */
	enum STRAND { FWD = 1, REV };

	/** default constructor */
	MEM() = default;

	/** construct an MEM with all info */
	MEM(const PrimarySeq* seq, const FMIndex* fmidx, STRAND strand,
			uint64_t from, uint64_t to, uint64_t SAstart, uint64_t SAend,
			const vector<Loc>& locs)
	: seq(seq), fmidx(fmidx), strand(strand), from(from), to(to), SAstart(SAstart), SAend(SAend), locs(locs)
	{  }

	/** construct an MEM with all info but not locs */
	MEM(const PrimarySeq* seq, const FMIndex* fmidx, STRAND strand,
			uint64_t from, uint64_t to, uint64_t SAstart, uint64_t SAend)
	: seq(seq), fmidx(fmidx), strand(strand), from(from), to(to), SAstart(SAstart), SAend(SAend)
	{  }

	/* member methods */
	/** get length of this MEM */
	uint64_t length() const {
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
		locs.clear();
	}

	/**
	 * fill the locs information of this MEM
	 * locs are lazy-filled to save computation time
	 */
	MEM& findLocs();

	/**
	 * get loglikelihood of this MEM
	 * @param  ignoreQual  whether to ignore quality even exists
	 * @return  log-pvalue of observing this MEM by random, using base-frequency and optionally base quality
	 */
	double loglik() const;

	/** get the pvalue of observing this MEM by random */
	double pvalue() const {
		return ::exp(loglik());
	}

	/** get the E-value of observing this MEM on a known size database */
	double evalue() const {
		return fmidx->length() * pvalue();
	}

	/* static member methods */
	/** test whether two MEM overlap on the seq */
	static bool isOverlap(const MEM& mem1, const MEM& mem2) {
		return mem1.from < mem2.to && mem1.to > mem2.from;
	}

	/** get the seq-distance of two mem,
	 * return 0 if they are overlapping
	 */
	static int64_t seqDist(const MEM& mem1, const MEM& mem2);

	/** get the DB-distance of two mem,
	 * return 0 if they are overlapping,
	 * or -1 both locs are empty
	 */
	static int64_t dbDist(const MEM& mem1, const MEM& mem2);

	/** get the number of indeals of two MEM
	 * return positive number if insertion, negative if deletion, or 0 if none
	 */
	static int64_t nindel(const MEM& mem1, const MEM& mem2) {
		return seqDist(mem1, mem2) - dbDist(mem2, mem2);
	}

	/** get the indel rate relative to the size of their mapped seq
	 * return positive number if insertion, negative if deletion, or 0 if none
	 */
	static double rindel(const MEM& mem1, const MEM& mem2) {
		return nindel(mem1, mem2) / static_cast<double> (mem1.seq->length());
	}

	/* member fields */
	const PrimarySeq* seq = nullptr;
	const FMIndex* fmidx = nullptr;
	STRAND strand = FWD;
	int64_t from = 0; /* 0-based relative start on seq */
	int64_t to = 0;   /* 1-based relative end on seq */
	int64_t SAstart = 0; /* 0-based start position on SA */
	int64_t SAend = 0;   /* 1-based end position on SA */
	vector<Loc> locs; /* all Loc this MEM matches to w/ reversed coordinates */

	/* static methods */
	/**
	 * get an MEM of a given seq starting at given position relative to the seq
	 * the locs will not be filled by this method
	 * @param seq  seq to search, must be in reversed orientation of this FM-index
	 * @param i  relative position of the seq
	 * @param strand  which strand to search
	 */
	static MEM findMEM(const PrimarySeq* seq, const FMIndex* fmidx, uint64_t from = 0, STRAND strand = FWD);
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_MEM_H_ */
