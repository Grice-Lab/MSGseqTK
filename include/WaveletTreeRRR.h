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
 *  A WaveletTreeRRR uses BitSeqRRR to store the bits internally for balanced speed and storage
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
	WaveletTreeRRR(const basic_string<uIntType>& src, size_t min = -1, size_t max = -1)
	: wid(sizeof(uIntType) * Wb), min(min), max(max) {
		if(this->min == -1)
			this->min = *std::min_element(src.begin(), src.end());
		if(this->max == -1)
			this->max = *std::max_element(src.begin(), src.end());

		build(src);
	}

	/* member methods */
	/** test whether the i-th bit of val is set */
	bool test(size_t val, size_t i) const {
		assert (i < height);
		return (val & (1UL << height - i - 1)) != 0;
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
	 * @return rank of s till position i
	 * @override  base class implementation
	 */
	virtual size_t rank(size_t s, size_t i) const;

	/**
	 * get the position of the r-th occurrence of symbol s
	 * @param s  symbol to search
	 * @param r  required occurrence/rank of s
	 * @return  position in this Seq, or -1 if i = 0, or len if i > total # of s
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

	/** get count of given symbol */
	virtual size_t getCount(size_t s) const {
		return OCC[s + 1] - OCC[s] + 1;
	}

	/** get cumulative count of given symbol */
	virtual size_t getOcc(size_t s) const {
		return OCC[s];
	}

	/**
	 * find the q-th smallest element in Seq[l..r]
	 */
	virtual size_t quantile(size_t left, size_t right, size_t q);

	/* find the q-th smallest element in T[l..r] and return the freq */
	virtual pair<size_t, size_t> quantile_freq(size_t left, size_t right, size_t q);

	/* utility methods */
private:
	/**
	 * build the basic fields
	 */
	template<typename uIntType>
	void build(const basic_string<uIntType>& src);

	/**
	 * recursively build Wavelet Tree bstrs of given level
	 * @param sym  (encoded) input symbol substring at given level
	 * @param level  level to build
	 * @param offset  offset relative to the original symbol
	 */
	template<typename uIntType>
	void build_level(vector<BitStr32>& bstrs, const basic_string<uIntType>& sym, size_t level, size_t offset = 0);

	/* member fields */
private:
	size_t height = 0; /* height of Wavelet Tree */
	size_t wid = 0; /* nBits of the original input string */
	size_t min = 0; /* min value of input symbols */
	size_t max = 0; /* max value of input symbols */
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
void WaveletTreeRRR::build(const basic_string<uIntType>& src) {
	/* build basic fields */
	n = src.length();
	sigma = max - min + 1;
	height = bits(max + 1);
	OCC.resize(max + 2); /* 1-based occurence is dummy position */
	/* get original OCC */
	OCC[0] = 0;
	for (uIntType ch : src)
        OCC[ch + 1]++; /* avoid zeros in OCC */

	/* enlarge symbols if required */
	size_t to_add = 0;
	for (size_t i = 1; i <= max + 1; ++i)
		if(OCC[i] == 0)
			to_add++;

	n += to_add;

	/* construct temporary bitstrgings */
	vector<BitStr32> bstrs;  /* an array of BitStr32 of length height, to store the intermediate status during construction */
//	bstrs.reserve(height);
	for(size_t i = 0; i < height; ++i)
		bstrs.push_back(BitStr32(n));

	if(to_add == 0)
		build_level(bstrs, src, 0);
	else {
		basic_string<uIntType> srcN = src; /* make a local copy */
		for (size_t i = 1; i <= max + 1; ++i) {
			if (OCC[i] == 0) {
				/* append a new character (i - 1) as OCC = 1 */
				OCC[i]++;
				srcN.push_back(i - 1);
			}
		}
		build_level(bstrs, srcN, 0);
	}
	cerr << "bitstrs built" << endl;

	/* build the BitSeqs from BitStrs */
//	bseqs.reserve(height);
	for(const BitStr32& bstr : bstrs) /* use a copy instead of const reference to prevent move construction */
		bseqs.push_back(BitSeqRRR(bstr));
	/* clear intermediate storates */
	cerr << "bitseqs built" << endl;

	/* build cumulative OCC */
	for(size_t i = 1; i <= max + 1; ++i)
		OCC[i] += OCC[i - 1];
}

template<typename uIntType>
void WaveletTreeRRR::build_level(vector<BitStr32>& bstrs, const basic_string<uIntType>& sym, size_t level, size_t offset) {
	if(level == height)
		return;

	size_t cleft = 0;
	for (size_t i = 0; i < n; ++i)
		if (!test(sym[i], level))
			cleft++;

	size_t cright = n - cleft;

	basic_string<uIntType> left;
	basic_string<uIntType> right;

	assert(bstrs[level].length() == n);

	left.reserve(cleft);
	right.reserve(cright);

	for (size_t i = 0; i < sym.length(); ++i) {
		bool flag = test(sym[i], level);
		bstrs[level].set(i + offset, flag);
		if(!flag)
			left.push_back(sym[i]);
		else
			right.push_back(sym[i]);
	}

	/* build level recursevely */
	build_level(bstrs, left, level + 1, offset);
	build_level(bstrs, right, level + 1, offset + cleft);
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* WAVELETTREERRR_H_ */
