/*
 * WaveletTreeRRR.h
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#ifndef WAVELETTREERRR_H_
#define WAVELETTREERRR_H_

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>
#include "Seq.h"
#include "BitStr.h"
#include "BitSeqRRR.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

using std::basic_string;
using std::vector;
using std::pair;

/**
 * A Raman, Raman and Rao's WaveletTree implementation of Seq, based on the algorithm described in
 *  [1] R. Raman, V. Raman and S. Rao. Succinct indexable dictionaries with applications
 *     to encoding $k$-ary trees and multisets. SODA02.
 *  [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 *  A WaveletTreeRRR uses BitSeqRRR to store the bits internally for balancing speed and storage
 */
class WaveletTreeRRR: public Seq {
public:
	/* constructors */
	/** default constructor */
	WaveletTreeRRR() = default;

	/**
	 * build a WaveletTree from a std::basic_string of any type in any alphabet
	 * @param src  copy of input string
	 */
	template<typename uIntType>
	WaveletTreeRRR(const basic_string<uIntType>& src, size_t min = -1, size_t max = -1, size_t sample_rate = BitSeqRRR::DEFAULT_SAMPLE_RATE)
	: wid(sizeof(uIntType) * Wb), min(min), max(max), sample_rate(sample_rate) {
		if(this->min == -1)
			this->min = *std::min_element(src.begin(), src.end());
		if(this->max == -1)
			this->max = *std::max_element(src.begin(), src.end());
		build(src);
	}

	/**
	 * build a WaveletTree from a std::basic_string of any type in any alphabet, non-const version
	 * @param src  copy of input string
	 */
	template<typename uIntType>
	WaveletTreeRRR(basic_string<uIntType>& src, size_t min = -1, size_t max = -1, size_t sample_rate = BitSeqRRR::DEFAULT_SAMPLE_RATE)
	: wid(sizeof(uIntType) * Wb), min(min), max(max), sample_rate(sample_rate) {
		if(this->min == -1)
			this->min = *std::min_element(src.begin(), src.end());
		if(this->max == -1)
			this->max = *std::max_element(src.begin(), src.end());
		build(src);
	}

	/* member methods */
	/** test whether the i-th bit of val is set */
	bool test(size_t val, size_t i) const {
		return val & (1UL << height - i - 1);
	}

	/**
	 * get the i-th symbol in this Seq
	 * @override the base class method
	 */
	virtual size_t access(size_t i) const;

	/**
	 * get the i-th bit of this BitSeq by testing the rank1 difference, and return the corresponding rank the same time
	 * @param i  pos
	 * @param r  rank of i if bit is on, or i - r + 1 if bit is off
	 * @return  symbol at pos i
	 * @override  base class method
	 */
	virtual size_t access(size_t i, size_t& r) const;

	/**
	 * get #occurrence of symbol s until position i (0-based, inclusive)
	 * @param s  symbol
	 * @param i  position
	 * @return rank of s till position i, or -1 if i >= n
	 * @override  base class implementation
	 */
	virtual size_t rank(size_t s, size_t i) const;

	/**
	 * get the position of the r-th occurrence of symbol s
	 * @param s  symbol to search
	 * @param r  required occurrence/rank of s
	 * @return  position in this Seq, or -1 if r = 0, or len if r > total # of s
	 * @override  base class implementation
	 */
	virtual size_t select(size_t s, size_t r) const;

	/**
	 * save this BitSeq to binary output
	 * @override  base class method
	 */
	virtual ostream& save(ostream& out) const;

	/**
	 * load data from a binary input
	 * @override  base class method
	 */
	virtual istream& load(istream& in);

	/** reset this Seq to default state
	 * @override  base class method
	 */
	virtual void reset() {
		bseqs.clear();
		OCC.clear();
		max = 0;
		min = 0;
		wid = 0;
		height = 0;
		Seq::reset();
	}

	/** get the size of the structure in bytes */
	virtual size_t getBytes() const;

	/** get count of given symbol */
	virtual size_t getCount(size_t s) const {
		return OCC[s + 1] - OCC[s] + 1;
	}

	/** get cumulative count of given symbol */
	virtual size_t getOcc(size_t s) const {
		return OCC[s];
	}

	/* non-member operators */
	friend bool operator==(const WaveletTreeRRR& lhs, const WaveletTreeRRR& rhs);

	/* utility methods */
private:
	/**
	 * build the basic fields from a copy of the input string
	 */
	template<typename uIntType>
	void build(const basic_string<uIntType>& src);

	/**
	 * build the basic fields from a copy of the input string, non-cost version
	 */
	template<typename uIntType>
	void build(basic_string<uIntType>& src);

	/**
	 * recursively build Wavelet Tree bstrs of given level
	 * @param sym  (encoded) input symbol substring at given level
	 * @param level  level to build
	 * @param offset  offset relative to the original symbol
	 */
	template<typename uIntType>
	void build_level(vector<BitStr32>& bstrs, const basic_string<uIntType>& sym,size_t level, size_t offset = 0);

	/**
	 * recursively build Wavelet Tree bstrs of given level, non-const version
	 * @param sym  (encoded) input symbol substring at given level
	 * @param level  level to build
	 * @param offset  offset relative to the original symbol
	 */
	template<typename uIntType>
	void build_level(vector<BitStr32>& bstrs, basic_string<uIntType>& sym,size_t level, size_t offset = 0);

	/* member fields */
private:
	size_t height = 0; /* height of Wavelet Tree */
	size_t wid = 0; /* nBits of the original input string */
	size_t min = 0; /* min value of input symbols */
	size_t max = 0; /* max value of input symbols */
	size_t sample_rate = 0; /* sample_rate for underlying BitSeqRRR structures */
	vector<size_t> OCC; /* 0-based cumulative occurence of each symbol in alphabet, with length max + 2 */
	vector<BitSeqRRR> bseqs; /* an array of BitSeqRRR, with length height */

	/* static methods */
public:
	static size_t get_start(size_t symbol, size_t mask) {
		return symbol & mask;
	}

	static size_t get_end(size_t symbol, size_t mask) {
		return get_start(symbol, mask) + ~mask + 1;
	}
};

template<typename uIntType>
inline void WaveletTreeRRR::build(const basic_string<uIntType>& src) {
	/* build basic fields */
	n = src.length();
	sigma = max - min + 1;
	height = bits(max);
	OCC.resize(max + 2); /* 1-based occurence is dummy position */

	/* get original OCC */
	OCC[0] = 0;
	for (uIntType ch : src)
        OCC[ch + 1]++; /* avoid zeros in OCC */

	/* construct intermediate BitStrs */
	vector<BitStr32> bstrs; /* intermediate BitStrs */
	bstrs.reserve(height);
	for(size_t i = 0; i < height; ++i)
		bstrs.push_back(BitStr32(n));

	/* build levels */
	build_level(bstrs, src, 0);

	/* build the BitSeqs from BitStrs */
	bseqs.reserve(height);
	for(const BitStr32& bs : bstrs)
		bseqs.push_back(BitSeqRRR(bs, sample_rate));

	/* build cumulative OCC */
	for(size_t i = 1; i <= max + 1; ++i)
		OCC[i] += OCC[i - 1];
}

template<typename uIntType>
inline void WaveletTreeRRR::build(basic_string<uIntType>& src) {
	/* build basic fields */
	n = src.length();
	sigma = max - min + 1;
	height = bits(max);
	OCC.resize(max + 2); /* 1-based occurence is dummy position */

	/* get original OCC */
	OCC[0] = 0;
	for (uIntType ch : src)
        OCC[ch + 1]++; /* avoid zeros in OCC */

	/* construct intermediate BitStrs */
	vector<BitStr32> bstrs; /* intermediate BitStrs */
	bstrs.reserve(height);
	for(size_t i = 0; i < height; ++i)
		bstrs.push_back(BitStr32(n));

	/* build levels */
	build_level(bstrs, src, 0);
	src.clear();
	src.shrink_to_fit();

	/* build the BitSeqs from BitStrs */
	bseqs.reserve(height);
	for(const BitStr32& bs : bstrs)
		bseqs.push_back(BitSeqRRR(bs, sample_rate));

	/* build cumulative OCC */
	for(size_t i = 1; i <= max + 1; ++i)
		OCC[i] += OCC[i - 1];
}

template<typename uIntType>
inline void WaveletTreeRRR::build_level(vector<BitStr32>& bstrs, const basic_string<uIntType>& sym,
		size_t level, size_t offset) {
	if(level == height || sym.empty())
		return;
	const size_t N = sym.length();
	BitStr32& bs = bstrs[level];

	assert(bs.length() == n);

	size_t nbits = 0;
	for (size_t i = 0; i < N; ++i) {
		if(test(sym[i], level)) {
			bs.set(i + offset);
			nbits++;
		}
	}

	basic_string<uIntType> left(N - nbits, 0);
	basic_string<uIntType> right(nbits, 0);
	for(size_t k = 0, i = 0, j = 0; k < N; ++k)
		if(test(sym[k], level))
			right[j++] = sym[k];
		else
			left[i++] = sym[k];

	/* build level recursevely */
	build_level(bstrs, left, level + 1, offset);
	build_level(bstrs, right, level + 1, offset + left.length());
}

template<typename uIntType>
inline void WaveletTreeRRR::build_level(vector<BitStr32>& bstrs, basic_string<uIntType>& sym,
		size_t level, size_t offset) {
	if(sym.empty())
		return;
	if(level == height) {
		sym.clear();
		sym.shrink_to_fit();
		return;
	}
	const size_t N = sym.length();
	BitStr32& bs = bstrs[level];

	assert(bs.length() == n);

	size_t nbits = 0;
	for (size_t i = 0; i < N; ++i) {
		if(test(sym[i], level)) {
			bs.set(i + offset);
			nbits++;
		}
	}

	const uIntType FLAG = std::numeric_limits<uIntType>::max();
//	basic_string<uIntType> left(N - nbits, 0);
	basic_string<uIntType> right(nbits, 0);
	for(size_t k = 0, j = 0; k < N; ++k) {
		if(test(sym[k], level)) {
			right[j++] = sym[k];
			sym[k] = FLAG;
		}
	}
	/* use sym itself as left partition */
	basic_string<uIntType>& left = sym;
	left.erase(std::remove(left.begin(), left.end(), FLAG), left.end());
	left.shrink_to_fit();

	/* build level recursevely */
	build_level(bstrs, left, level + 1, offset);
	build_level(bstrs, right, level + 1, offset + N - nbits);
}

inline bool operator==(const WaveletTreeRRR& lhs, const WaveletTreeRRR& rhs) {
	return dynamic_cast<const Seq&>(lhs) == dynamic_cast<const Seq&>(rhs) &&
			lhs.height == rhs.height && lhs.wid == rhs.wid &&
			lhs.min == rhs.min && lhs.max == rhs.max &&
			lhs.OCC == rhs.OCC && lhs.bseqs == rhs.bseqs;
}

inline bool operator!=(const WaveletTreeRRR& lhs, const WaveletTreeRRR& rhs) {
	return !(lhs == rhs);
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* WAVELETTREERRR_H_ */
