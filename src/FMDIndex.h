/*
 * FMDIndex.h
 *  A bidirectional FM-index for both the forward and reverse-complement of a DNAseq
 *  proposed in Li Heng Bioinformatics 2012
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#ifndef SRC_FMDINDEX_H_
#define SRC_FMDINDEX_H_

#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include "DNAseq.h"
#include "QualStr.h"
#include "GLoc.h"
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
using EGriceLab::libSDS::Seq;
using EGriceLab::libSDS::WaveletTreeRRR;

/**
 * FMD-index with bi-directinal search/extend methods
 */
class FMDIndex {
public:
	/* typedefs */
	typedef array<int64_t, DNAalphabet::SIZE> BCarray_t; /* fixed array to store base counts */
	typedef vector<int64_t> SArray_t; /* store sampled Suffix-Array in std::vector */

	/* constructors */
	/** Default constructor */
	FMDIndex() = default;

	/** Construct an FMIndex from a given DNAseq */
	explicit FMDIndex(const DNAseq& seq, bool keepSA = false)
	: keepSA(keepSA) {
		build(seq);
	}

	/** construct an FMIndex from pre-built values */
	FMDIndex(const BCarray_t& B, const BCarray_t& C, const DNAseq& bwtSeq, bool keepSA = false);

	/**
	 * LF-mapping given position and character
	 * @param i  0-based location in bwt
	 * @param b  base/character for lookup
	 * @return  1-based index on F column (original seq)
	 */
	int64_t LF(uint8_t b, int64_t i) const {
		return C[b] + bwt.rank(b, i);
	}

	/**
	 * LF-mapping at given position
	 * @param i  0-based location in bwt
	 * @return 1-based index on F column (original seq)
	 */
	int64_t LF(int64_t i) const {
		return LF(bwt.access(i), i);
	}

	/** test whether this RRFMIndex is initiated */
	bool isInitiated() const {
		return !bwt.empty();
	}

	/** check whether this FMDIndex is bi-directional by checking counts of complement bases */
	bool isBiDirectional() const;

	/** get the length of this index */
	int64_t length() const {
		return bwt.length();
	}

	/** get baseCount array of this index */
	const BCarray_t& getBaseCount() const {
		return B;
	}

	/** get baseCount of given base in fwd seq */
	int64_t getBaseCount(uint8_t b) const {
		return B[b];
	}

	/** get baseCount of basic bases (A,T,C,G) */
	int64_t getBasicBaseCount() const {
		return B[DNAalphabet::A] + B[DNAalphabet::C] + B[DNAalphabet::G] + B[DNAalphabet::T];
	}

	/** get baseCount of extended/ambigous bases (IUPAC non-A,T,C,G and non-gap) */
	int64_t getExtBaseCount() const {
		return getCumCount(DNAalphabet::NT16_MAX) - getBasicBaseCount() - B[0];
	}

	/** get cumulative base count of given base */
	int64_t getCumCount(uint8_t b) const {
		return C[b];
	}

	/** get total bases in this FMD-index, alias to length() */
	const int64_t totalBases() const {
		return length();
	}

	/**
	 * get the frequency of a given base
	 * @return  base frequency if is a basic base, or the mapped basic base frequency if is an IMPAC extension
	 */
	double getBaseFreq(uint8_t b) const {
		return static_cast<double>(getBaseCount(DNAalphabet::toBasic(b))) / length();
	}

	/** get loglik() of a given base as the log base-frequency */
	double loglik(uint8_t b) const {
		return std::log(getBaseCount(DNAalphabet::toBasic(b))) - std::log(length());
	}

	/**
	 * Build an FMDIndex from a combined seq, in which seq is always in the order of R0R0'R1R1', etc
	 * @param seq pre-combined bidirectional seq, must be null (GAP_BASE) terminated
	 * @return a fresh allocated FMDIndex
	 * @throw std::length_error if seq length is too long
	 */
	FMDIndex& build(const DNAseq& seq);

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
	int64_t accessSA(int64_t i) const;

	/** build SAbit and SAsampled from the internal BWT, where raw BWTseq not available */
	void buildSA();

protected:
	/** build SAbit and SAsampled from a given raw BWTseq, which is available during merge operations */
	void buildSA(const DNAseq& bwtSeq);

	/**
	 * build SAbit and SAsampled from given SA and raw BWTseq, which are available during direct construction
	 * this version of buildSA is parallelized by openMP
	 */
	void buildSA(const int64_t* SA, const DNAseq& bwtSeq);

	/** build interleaving BitVector for two FMD-index, use parallelization optionally */
	static BitStr32 buildInterleavingBS(const FMDIndex& fmdidx1, const FMDIndex& fmdidx2);

	/** merge two DNAseq by an interleaving bitvector, use parallelization optionally */
	static DNAseq mergeBWT(const WaveletTreeRRR& bwt1, const WaveletTreeRRR& bwt2, const BitStr32& bstr);

public:
	/**
	 * backward extension of a bi-interval [p, q, s]
	 * @return  the new size as an indication of success or not
	 */
	int64_t backExt(int64_t& p, int64_t& q, int64_t& s, uint8_t b) const;

	/**
	 * forward extension of a bi-interval [p, q, s]
	 * @return  the new size as an indication of success or not
	 */
	int64_t fwdExt(int64_t& p, int64_t& q, int64_t& s, uint8_t b) const {
		return backExt(q, p, s, DNAalphabet::complement(b));
	}

	/**
	 * count present times of a DNAseq pattern in forward version
	 */
	int64_t count(const DNAseq& pattern) const;

	/** merge this RRFMIndex with another index with only very little overhead memory
	 * using the BWT-merge algorithm described in
	 * Burrows-Wheeler transform for terabases, Jouni Sirén, 2016 Data Compression Conference
	 */
	FMDIndex& operator+=(const FMDIndex& other);

	/** get the encoded BWT of the original seq */
	DNAseq getBWT() const;

	/** get the BWT string of the original seq */
	string getBWTStr() const {
		return dna::decode(getBWT());
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
	vector<GLoc> locateAll(const DNAseq& pattern, GLoc::STRAND strand = GLoc::FWD) const {
		return strand == GLoc::FWD ? locateAllFwd(pattern) : locateAllRev(pattern);
	}

	/** locate all matches in FWD direction */
	vector<GLoc> locateAllFwd(const DNAseq& pattern) const;

	vector<GLoc> locateAllRev(const DNAseq& pattern) const;

	/**
	 * reverse a loc on this FM-index
	 * @param i  0-based loc
	 * @return  1-based loc on the reversed direction
	 */
	int64_t reverseLoc(int64_t i) const {
		return length() - 1 - i;
	}

	/** reverse the coordinates of a Loc */
	Loc reverseLoc(const Loc& loc) const {
		return Loc(length() - 1 - loc.getEnd(), length() - 1 - loc.getStart()); // FM-index length is 1 more than the original seq length
	}

	/* non-member functions */
	friend FMDIndex operator+(const FMDIndex& lhs, const FMDIndex& rhs);

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
	BCarray_t B = { };  // combined base count
	BCarray_t C = { };  // combined cumulative count
	WaveletTreeRRR bwt; /* Wavelet-Tree transformed BWT string for combined seq */
	bool keepSA = false; /* whether to keep SAidx and SAsampled during building */
	BitSeqRRR SAbit; /* BitSeq index telling whether a SA was sampled */
	SArray_t SAsampled; /* sampled SA vector */

	/* static member fields */
public:
	static const int64_t MAX_LENGTH = std::numeric_limits<int64_t>::max();
};

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_FMDINDEX_H_ */
