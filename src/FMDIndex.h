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
using std::basic_string;
using std::array;
using EGriceLab::libSDS::BitStr32;
using EGriceLab::libSDS::BitSeqRRR;
using EGriceLab::libSDS::Seq;
using EGriceLab::libSDS::WaveletTreeRRR;

/**
 * FMD-index with bi-directinal search/extend methods
 * a fully built FMD-index keeps the base counts, a WaveletTreeRRR compressed BWT,
 * a BitSeqRRR compressed SA index and sampling,
 * and optionall the uncompressed BWT, which is useful for faster SA build algorithm during database construction
 */
class FMDIndex {
public:
	/* typedefs */
	typedef array<int64_t, DNAalphabet::SIZE> BCarray_t; /* fixed array to store base counts */
//	typedef vector<int64_t> SArray_t; /* store sampled Suffix-Array in std::vector */
	typedef basic_string<int64_t> SAarr_t; /* store sampled Suffix-Array values */

	/* constructors */
	/** Default constructor */
	FMDIndex() = default;

	/** Construct an FMDIndex from a given DNAseq
	 * it will build the counts and BWT,
	 * and optionally build the SA (SAidx and SAsampled)
	 */
	explicit FMDIndex(const DNAseq& seq, bool keepSA = true);

	/* member methods */
	/** test whether contains uncompressed BWT */
	bool hasBWT() const {
		return !bwt.empty();
	}

	/** test whether contains SA */
	bool hasSA() const {
		return !SAidx.empty();
	}

	/** get the encoded BWT of the original seq */
	DNAseq getBWT() const;

	/** clear uncompressed BWT */
	FMDIndex& clearBWT() {
		bwt.clear();
		return *this;
	}

	/** clear SA */
	FMDIndex& clearSA() {
		SAidx.reset();
		SAsampled.clear();
//		SAgap.clear();
		return *this;
	}

	/** get the original seq */
	DNAseq getSeq() const;

	/** access BWT at given position */
	nt16_t accessBWT(int64_t i) const {
		return hasBWT() ? bwt[i] : bwtRRR.access(i);
	}

	/*
	 * Access SA at given position, either by directly searching the stored value or the next sampled value
	 * @param i  0-based on SA
	 * @return  0-based position on original seq
	 */
	int64_t accessSA(int64_t i) const;

	/**
	 * LF-mapping given position and character
	 * @param i  0-based location in bwt
	 * @param b  base/character for lookup
	 * @return  1-based index on F column (original seq)
	 */
	int64_t LF(nt16_t b, int64_t i) const {
		return C[b] + bwtRRR.rank(b, i);
	}

	/**
	 * LF-mapping at given position
	 * @param i  0-based location in bwt
	 * @return 1-based index on F column (original seq)
	 */
	int64_t LF(int64_t i) const {
		return LF(accessBWT(i), i);
	}

	/** test whether this RRFMIndex is initiated */
	bool isInitiated() const {
		return !bwtRRR.empty();
	}

	/** check whether this FMDIndex is bi-directional by checking counts of complement bases */
	bool isBiDirectional() const;

	/** get the length of this index */
	int64_t length() const {
		return bwtRRR.length();
	}

	/** get baseCount array of this index */
	const BCarray_t& getBaseCount() const {
		return B;
	}

	/** get baseCount of given base in fwd seq */
	int64_t getBaseCount(nt16_t b) const {
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
	int64_t getCumCount(nt16_t b) const {
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
	double getBaseFreq(nt16_t b) const {
		return static_cast<double>(getBaseCount(DNAalphabet::toBasic(b))) / length();
	}

	/** get loglik() of a given base as the log base-frequency */
	double loglik(nt16_t b) const {
		return std::log(getBaseCount(DNAalphabet::toBasic(b))) - std::log(length());
	}

	/**
	 * save raw object data to output
	 */
	ostream& save(ostream& out) const;

	/**
	 * load raw object data from input
	 */
	istream& load(istream& in);

	/** build SAidx and SAsampled */
	FMDIndex& buildSA();

public:
	/**
	 * backward extension of a bi-interval [p, q, s]
	 * @return  the new size as an indication of success or not
	 */
	int64_t backExt(int64_t& p, int64_t& q, int64_t& s, nt16_t b) const;

	/**
	 * forward extension of a bi-interval [p, q, s]
	 * @return  the new size as an indication of success or not
	 */
	int64_t fwdExt(int64_t& p, int64_t& q, int64_t& s, nt16_t b) const {
		return backExt(q, p, s, DNAalphabet::complement(b));
	}

	/**
	 * count present times of a DNAseq pattern in forward version
	 */
	int64_t count(const DNAseq& pattern) const;

	/** merge this RRFMIndex with another index with only very little overhead memory
	 * using the BWT-merge algorithm described in
	 * Burrows-Wheeler transform for terabases, Jouni Sirén, 2016 Data Compression Conference
	 * it does not update the SA (SAidx and SAsampled) automatically to save unnecessary update,
	 * you will need to call buildSA() yourself
	 */
	FMDIndex& operator+=(const FMDIndex& other);

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

protected:
	/* helper methods */
	/**
	 * build alphabet counts from a DNAseq
	 */
	FMDIndex& buildCounts(const DNAseq& seq);

	/** build BWT uncompressed and compressed,
	 * and optionally build the SA */
	FMDIndex& buildBWT(const DNAseq& seq, bool keepSA = true);

	/**
	 * build SAidx and SAsampled with given internal SA, which are available during direct construction
	 */
	FMDIndex& buildSA(const int64_t* SA);

	/** build interleaving BitVector for two FMD-index, use parallelization optionally */
	static BitStr32 buildInterleavingBS(const FMDIndex& lhs, const FMDIndex& rhs);

	/** merge two DNAseq by an interleaving bitvector, use parallelization optionally */
	static DNAseq mergeBWT(const FMDIndex& lhs, const FMDIndex& rhs, const BitStr32& bstrM);

	/* member fields */
private:
	BCarray_t B = { };  // combined base count
	BCarray_t C = { };  // combined cumulative count
	DNAseq bwt; // uncompressed bwt
	WaveletTreeRRR bwtRRR; /* Wavelet-Tree transformed BWT string for combined seq */
	BitSeqRRR SAidx; /* BitSeq index telling whether a SA was sampled */
	SAarr_t SAsampled; /* sampled SA */

	/* static member fields */
public:
	static const int RRR_SAMPLE_RATE = 32; /* RRR sample rate for BWT */
	static const int SA_SAMPLE_RATE = 32;  /* sample rate for SA */
	static const int64_t MAX_LENGTH = std::numeric_limits<int64_t>::max();
};

inline FMDIndex operator+(const FMDIndex& lhs, const FMDIndex& rhs) {
	FMDIndex fmdidxM(lhs);
	fmdidxM += rhs;
	return fmdidxM;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_FMDINDEX_H_ */
