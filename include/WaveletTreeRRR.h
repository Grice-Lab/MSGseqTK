/*
 * WaveletTreeRRR.h
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#ifndef WAVELETTREERRR_H_
#define WAVELETTREERRR_H_

#include <string>
#include <algorithm>
#include "Seq.h"
#include "BitStr.h"
#include "BitSeqRRR.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

using std::basic_string;
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

	/** destructor */
	virtual ~WaveletTreeRRR() {
		delete[] OCC;
		for(unsigned int i = 0; i < height; ++i)
			delete[] bseqs[i];
		delete[] bseqs;
	}

	/** copy constructor */
	WaveletTreeRRR(const WaveletTreeRRR& other);

	/** copy assignment operator using copy-swap */
	WaveletTreeRRR& operator=(WaveletTreeRRR other) {
		return swap(other);
	}

	/**
	 * build a WaveletTree from a std::basic_string of any type in any alphabet
	 * @param src  copy of input string
	 */
	template<typename uIntType>
	WaveletTreeRRR(const basic_string<uIntType>& src, unsigned int min = -1, unsigned int max = -1) : min(min), max(max) {
		const size_t wSrc = sizeof(uIntType) * Wb;
		if(this->min == -1)
			this->min = *std::min_element(src.begin(), src.end());
		if(this->min == -1)
			this->min = *std::max_element(src.begin(), src.end());
		height = bits(this->max - this->min + 1);
		OCC = new size_t[sigma + 1](); /* Zero-initiation, 0 is dummy position */
		for (size_t i = 0; i < n; ++i)
            OCC[src[i] + 1]++; /* avoid zeros in OCC */

		BitStr32 sym(src); /* construct a BitStr copy of src */

		/* check how many we need to enlarge symbols */
		unsigned int to_add = 0;
		for (unsigned int i = 1; i <= max + 1; ++i)
			if(OCC[i] == 0)
				to_add++;

		sym.resize((src.length() + to_add) * wSrc); // resize sym, may not need at all

		n = src.length() + to_add;
		sigma = this->max - this->min + 1;

		if(to_add == 0) /* no need to modify */
			build_basic(src);
		else {
			basic_string<uIntType> sym(src); /* make a copy */
			for (unsigned int i = 1; i <= max + 1; ++i) {
				if (OCC[i] == 0) {
					/* append a new character (i - 1) as OCC = 1 */
					OCC[i]++;
					sym.push_back(i - 1);
				}
			}
			build_basic(sym);
		}
	}

	/* member methods */
	/** test whether the i-th bit of val is set */
	bool test(size_t val, unsigned int i) const {
		assert (i < height);
		return (val & (1 << height - i - 1)) != 0;
	}

	/**
	 * get the i-th symbol in this Seq
	 * @override the base class method
	 */
	virtual unsigned int access(size_t i) const;

	/**
	 * get #occurrence of symbol s until position i (0-based, inclusive)
	 * @param s  symbol
	 * @param i  position
	 * @return rank of s till position i
	 * @override  base class implementation
	 */
	virtual size_t rank(unsigned int s, size_t i) const;

	/**
	 * get the position of the r-th occurrence of symbol s
	 * @param s  symbol to search
	 * @param r  required occurrence/rank of s
	 * @return  position in this Seq, or -1 if i = 0, or len if i > total # of s
	 * @override  base class implementation
	 */
	virtual size_t select(unsigned int s, size_t r) const;

	/** copy this object with another */
	WaveletTreeRRR& swap(WaveletTreeRRR& other);

	/* utility methods */
private:
	/**
	 * build the basic fields from BitStr32 representing the original string
	 * @param bstr  a BitStr32 copy of input string
	 */
	void build_basic(const BitStr32& bstr);

	/**
	 * recursively build Wavelet Tree bstrs of given level
	 * @param bstr  a BitStr32 copy of the input string at given level
	 * @param level  level to build
	 *
	 */
	void build_level(const BitStr32& sym, unsigned int level = 0);

	/* member fields */
private:
	unsigned int height = 0; /* height of Wavelet Tree */
	unsigned int min = 0; /* min value of input symbols */
	unsigned int max = 0; /* max value of input symbols */
	size_t* OCC = nullptr; /* 0-based occurence of each symbol in alphabet, with length max + 2 */
	BitSeqRRR* bseqs = nullptr; /* an array of BitStrings, with length height */
};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* WAVELETTREERRR_H_ */
