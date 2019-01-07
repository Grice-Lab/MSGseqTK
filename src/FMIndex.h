/*
 * RFMIndex.h
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#ifndef SRC_FMINDEX_H_
#define SRC_FMINDEX_H_

#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include "DNAseq.h"
#include "QualStr.h"
#include "Loc.h"
#include "divsufsort_private.h"
#include "BitStr.h"
#include "BitSeqRRR.h"
#include "WaveletTreeRRR.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::vector;
using std::array;
using EGriceLab::libSDS::BitStr32;
using EGriceLab::libSDS::BitSeqRRR;
using EGriceLab::libSDS::WaveletTreeRRR;

/**
 * FM-index with merge support
 */
class FMIndex {
	typedef array<saidx_t, DNAalphabet::SIZE> BCarray_t; /* fixed array to store base counts */
	typedef vector<saidx_t> SArray_t; /* store sampled Suffix-Array in std::vector */

public:
	/* constructors */
	/** Default constructor */
	FMIndex() = default;

	/** Construct an FMIndex from a given (large) DNAseq */
	explicit FMIndex(const DNAseq& seq, bool keepSA = false)
	: keepSA(keepSA) {
		build(seq);
	}

	/** construct an FMIndex from pre-built values */
	FMIndex(const BCarray_t& B, const BCarray_t& C, const DNAseq& bwtSeq, bool keepSA = false);

	/**
	 * LF-mapping given position and character
	 * @param i  0-based location in bwt
	 * @param b  base/character for lookup
	 * @return  1-based index on F column (original seq)
	 */
	saidx_t LF(sauchar_t b, saidx_t i) const {
		return C[b] + bwt.rank(b, i);
	}

	/**
	 * LF-mapping at given position
	 * @param i  0-based location in bwt
	 * @return 1-based index on F column (original seq)
	 */
	saidx_t LF(saidx_t i) const {
		return LF(bwt.access(i), i);
	}

	/** test whether this RRFMIndex is initiated */
	bool isInitiated() const {
		return !bwt.empty();
	}

	/** get the length of this index */
	saidx_t length() const {
		return bwt.length();
	}

	/** get baseCount array of this index */
	const BCarray_t& getBaseCount() const {
		return B;
	}

	/** get baseCount of given base */
	saidx_t getBaseCount(sauchar_t b) const {
		return B[b];
	}

	/** get baseCount of basic bases (A,T,C,G) */
	saidx_t getBasicBaseCount() const {
		return B[DNAalphabet::A] + B[DNAalphabet::C] + B[DNAalphabet::G] + B[DNAalphabet::T];
	}

	/** get baseCount of extended/ambigous bases (IUPAC non-A,T,C,G) */
	saidx_t getExtBaseCount() const {
		return C[DNAalphabet::NT16_MAX] - getBasicBaseCount();
	}

	/** get cumulative base count of given base */
	saidx_t getCumCount(sauchar_t b) const {
		return C[b];
	}

	/** get total bases in this FM-index excluding gaps */
	const saidx_t totalBases() const;

	/**
	 * Build an RRFMIndex from a MSA object, old data is removed
	 * @param msa  pointer to an MSA object
	 * @return a fresh allocated RRFMIndex
	 * @throw std::length_error if seq length is too long
	 */
	FMIndex& build(const DNAseq& seq);

	/**
	 * save raw object data to output
	 */
	ostream& save(ostream& out) const;

	/**
	 * load raw object data from input
	 */
	istream& load(istream& in);

	/*
	 * Access SA at given position, either by directly searching the stored value or the next sampled value
	 * @param i  0-based on SA
	 * @return  0-based position on original seq
	 */
	saidx_t accessSA(saidx_t i) const;

	/** build SAbit and SAsampled from the internal BWT, where raw BWTseq not available */
	void buildSA();

protected:
	/** build SAbit and SAsampled from a given raw BWTSeq, which is available during many operations */
	void buildSA(const DNAseq& bwtSeq);

public:
	/**
	 * count present times of a DNAseq pattern in forward version
	 */
	saidx_t count(const DNAseq& pattern) const;

	/** merge this RRFMIndex with another index with only very little overhead memory
	 * using the BWT-merge algorithm described in
	 * Burrows-Wheeler transform for terabases, Jouni Sirén, 2016 Data Compression Conference
	 * Note that the other will be used as seq 1, to improve the performance of the algorithm
	 */
	FMIndex& operator+=(const FMIndex& other);

	/** get the encoded BWT of the original seq */
	DNAseq getBWT() const;

	/** get the BWT string of the original seq */
	string getBWTStr() const {
		return getBWT().toString();
	}

	/** get the original seq */
	DNAseq getSeq() const;

	bool hasSA() const {
		return keepSA;
	}

	/** clear stored SA information, if any */
	void clearSA() {
		SAsampled.clear();
		SAbit.reset();
		keepSA = false;
	}

	/** locate all matches to given pattern */
	vector<Loc> locateAll(const DNAseq& pattern) const;

	/**
	 * reverse a loc on this FM-index
	 * @param i  0-based loc
	 * @return  1-based loc on the reversed direction
	 */
	saidx_t reverseLoc(saidx_t i) const {
		return length() - 1 - i;
	}

	/** reverse the coordinates of a Loc on this FM-index */
	Loc reverseLoc(const Loc& loc) const {
		return Loc(length() - 1 - loc.end, length() - 1 - loc.start); // FM-index length is 1 more than the original seq length
	}

	/* non-member functions */
	friend FMIndex operator+(const FMIndex& lhs, const FMIndex& rhs);

private:
	/**
	 * build alphabet counts from a DNAseq
	 */
	void buildCounts(const DNAseq& seq);

	/** build BWT on the reversed string of seq */
	void buildBWT(const DNAseq& seq);

	static const int RRR_SAMPLE_RATE = 32; /* RRR sample rate for BWT */
	static const int SA_SAMPLE_RATE = 32;  /* sample rate for SA */

	/* member fields */
private:
	BCarray_t B = { }; /* base count of each alphabetbase */
	BCarray_t C = { }; /* cumulative count of each alphabet base, C[i] = B(0) + B(1) + ... + B(i-1) */
	WaveletTreeRRR bwt; /* Wavelet-Tree transformed BWT string for reversed DNAseq */
	bool keepSA = false; /* whether to keep SAidx and SAsampled during building */
	BitSeqRRR SAbit; /* BitSeq index telling whether a SA was sampled */
	SArray_t SAsampled; /* sampled SA vector */

	/* static member fields */
public:
	static const saidx_t MAX_LENGTH = std::numeric_limits<saidx_t>::max();
};

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_FMINDEX_H_ */
