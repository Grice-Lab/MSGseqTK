/*
 * RFMIndex.h
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#ifndef SRC_RFMINDEX_H_
#define SRC_RFMINDEX_H_

#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstdint> // C++11
#include <limits>
#include <stdexcept>
#include <memory>
#include "DNAseq.h"
#include "divsufsort.h" /* use 32bit version */
#include "WaveletTreeNoptrs.h"
#include "BitSequence.h"

namespace EGriceLab {
namespace MSGseqClean {

/**
 * A reduced FM-index with no locate function
 * but with merging support
 */
class RFMIndex {
	typedef cds_static::WaveletTreeNoptrs BWTRRR;
	typedef std::shared_ptr<BWTRRR> BWTRRR_ptr;

public:
	/* constructors */
	/** Default constructor */
	RFMIndex() = default;

	/** Construct a RRFMIndex from a given (large) DNAseq */
	explicit RFMIndex(const DNAseq& seq) {
		build(seq);
	}

	/** destructor */
	virtual ~RFMIndex() { 	}

	/**
	 * LF-mapping given position and character
	 * @param i  0-based location in bwt
	 * @param b  base/character for lookup
	 * @return  1-based index on F column (original seq)
	 */
	saidx_t LF(sauchar_t b, saidx_t i) const {
		return C[b] + bwt->rank(b, i);
	}

	/**
	 * LF-mapping at given position
	 * @param i  0-based location in bwt
	 * @return 1-based index on F column (original seq)
	 */
	saidx_t LF(saidx_t i) const {
		return LF(bwt->access(i), i);
	}

	/** test whether this RRFMIndex is initiated */
	bool isInitiated() const {
		return bwt != nullptr;
	}

	/** get the length of this index */
	saidx_t length() const {
		return bwt->getLength();
	}

	/**
	 * Build an RRFMIndex from a MSA object, old data is removed
	 * @param msa  pointer to an MSA object
	 * @return a fresh allocated RRFMIndex
	 * @throw std::length_error if seq length is too long
	 */
	RFMIndex& build(const DNAseq& seq);

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

	/** merge this RRFMIndex with another index with only very little overhead memory
	 * using the BWT-merge algorithm described in
	 * Burrows-Wheeler transform for terabases, Jouni Sirén, 2016 Data Compression Conference
	 */
	RFMIndex& operator+=(const RFMIndex& other);

	/** get the BWT string of the original seq */
	string getBWT() const;

	/** get the original seq */
	DNAseq getSeq() const;

	/* non-member functions */
	friend RFMIndex operator+(const RFMIndex& lhs, const RFMIndex& rhs);

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
	saidx_t C[UINT8_MAX + 1] = { 0 }; /* cumulative count of each alphabet frequency, with C[0] as dummy position */
	// DNAseq seq (concatenated) DNAseq underlying this index
	// uint32_t* SA = nullptr; /* 1-based sampled SA of reversed DNAseq */
	// cds_static::BitSequence* SAidx = nullptr; /* 1-based bit index for telling whether this SA position is sampled */
	BWTRRR_ptr bwt; /* Wavelet-Tree transformed BWT string for reversed DNAseq */

public:
	static const saidx_t MAX_LENGTH = std::numeric_limits<saidx_t>::max();
	static const char TERMINAL_SYMBOL = '$'; /* traditional terminal symbol of BWT */
};

inline RFMIndex operator+(const RFMIndex& lhs, const RFMIndex& rhs) {
	RFMIndex rrfm(lhs);
	return rrfm += rhs;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_RFMINDEX_H_ */
