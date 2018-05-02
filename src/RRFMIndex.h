/*
 * RRFMIndex.h
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#ifndef SRC_RRFMINDEX_H_
#define SRC_RRFMINDEX_H_

#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstdint> // C++11
#include <limits>
#include <stdexcept>
#include "DNAseq.h"
#include "divsufsort.h" /* use 32bit version */
#include "WaveletTreeNoptrs.h"
#include "BitSequence.h"

namespace EGriceLab {
namespace MSGseqClean {

/**
 * A Reversed and Reduced FM-index for ultra-fast seed query of large genome sequences
 */
class RRFMIndex {
public:
	/* constructors */
	/** Default constructor */
	RRFMIndex() = default;

	/** Construct a RRFMIndex from a given (large) DNAseq */
	explicit RRFMIndex(const DNAseq& seq) throw(std::length_error) {
		build(seq);
	}

	/* disable copy and assignment constructors */
	RRFMIndex(const RRFMIndex&) = delete;
	RRFMIndex& operator=(const RRFMIndex&) = delete;

	/** distructor */
	virtual ~RRFMIndex() {
		clear();
	}

	/** clear old data */
	virtual void clear() {
		delete bwt;
	}

	/**
	 * Build an RRFMIndex from a MSA object, old data is removed
	 * @param msa  pointer to an MSA object
	 * @return a fresh allocated RRFMIndex
	 * @throw std::length_error if seq length is too long
	 */
	RRFMIndex& build(const DNAseq& seq) throw(std::length_error);

	/**
	 * save raw object data to output
	 */
	ostream& save(ostream& out) const;

	/**
	 * load raw object data from input
	 */
	istream& load(istream& in);

	/**
	 * count present times of a DNAseq pattern in forward version
	 */
	saidx_t count(const DNAseq& pattern) const;

private:
	/**
	 * build alphabet counts from a DNAseq
	 */
	void buildCounts(const DNAseq& seq);

	/** build BWT on the reversed string of seq */
	void buildBWT(const DNAseq& seq);

	static const unsigned RRR_SAMPLE_RATE = 8; /* RRR sample rate for BWT */

	/* member fields */
private:
	saidx_t C[UINT8_MAX + 1]; /* cumulative count of each alphabet frequency, with C[0] as dummy position */
	// DNAseq seq (concatenated) DNAseq underlying this index
	// uint32_t* SA = nullptr; /* 1-based sampled SA of reversed DNAseq */
	// cds_static::BitSequence* SAidx = nullptr; /* 1-based bit index for telling whether this SA position is sampled */
	cds_static::WaveletTreeNoptrs* bwt = nullptr; /* Wavelet-Tree transformed BWT string for reversed DNAseq */

public:
	static const saidx_t MAX_LENGTH = std::numeric_limits<saidx_t>::max();

};

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_RRFMINDEX_H_ */
