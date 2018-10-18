/*
 * WaveletTreeRRR.h
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#ifndef WAVELETTREERRR_H_
#define WAVELETTREERRR_H_

#include "Seq.h"
#include "Mapper.h"
#include "BitStr.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

/**
 * A Raman, Raman and Rao's WaveletTree implementation of Seq, based on the algorithm described in
 *  [1] R. Raman, V. Raman and S. Rao. Succinct indexable dictionaries with applications
 *     to encoding $k$-ary trees and multisets. SODA02.
 *  [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 */
class WaveletTreeRRR: public Seq {
public:
	typedef size_t size_type;
	typedef BitStr<uint32_t> BitStrUint; /* default data type */
	/* constructors */
	/** default constructor */
	WaveletTreeRRR() = default;

	/** destructor */
	virtual ~WaveletTreeRRR() {
		delete[] OCC;
	}

	/** copy constructor */
	WaveletTreeRRR(const WaveletTreeRRR& other);

	/**
	 * build a WaveletTree from a C-string, the alphabet of src is assumbed to be compact with no zero-count alphabets
	 * @param src  input C-str
	 * @param n  length of input C-str
	 */
	template<typename uIntType>
	WaveletTreeRRR(const uIntType* src, size_type n, unsigned int maxV = 0) : n(n), max(maxV) {
		if(max == 0)
			max = max_value(src, n);
		height = bits(max);
		OCC = new size_t[max + 2](); /* Zero-initiation */
		for (size_t i = 0; i < n; ++i)
            OCC[src[i] + 1]++; /* avoid zeros in OCC */

		BitStr<uint32_t> bstr(src, n);

//		uint* new_symb = new uint[n + to_add];
        uint * new_symb = new uint[((n+to_add)*width)/W + 1]; // fixed by QZ

		/* copy old symbo to new symbols array */
        for(uint i = 0; i < n; ++i)
            set_field(new_symb, width, i, get_field(symbols, width, i)); // fixed by QZ
/*		for (uint i = 0; i < n; i++)
			new_symb[i] = symbols[i];*/

		if (deleteSymbols) {
			delete [] symbols;
			symbols = 0;
		}

		to_add = 0;
		for (uint i = 1; i <= max_v + 1; i++)
		if (OCC[i] == 0) {
			OCC[i]++;
            set_field(new_symb, width, n+to_add, i - 1); // fixed by QZ
			//new_symb[n + to_add] = i - 1;
			to_add++;
		}

		uint new_n = n + to_add;
		for(uint i = 1; i <= max_v + 1; i++)
			OCC[i] += OCC[i - 1];
		this->n = new_n;

		uint **_bm = new uint*[height];
		for(uint i = 0; i < height; i++) {
			_bm[i] = new uint[new_n / W + 1](); // zero-initiation fixed by QZ
/*			for(uint j = 0; j < new_n / W + 1; j++)
				_bm[i][j] = 0;*/
		}

		build_level(_bm, new_symb, width, 0, new_n, 0);
		bitstring = new BitSequence*[height];
		for(uint i=0;i< height; i++) {
			bitstring[i] = bmb->build(_bm[i], new_n);
			delete [] _bm[i];
		}
		delete [] _bm;

		if (!deleteSymbols)
			for (uint i = 0; i < n; i++)
				set_field(symbols, width, i, am->unmap(get_field(symbols, width, i)));
	}

	/** build a WaveletTree from a C-string and Mapper */
	template<typename uIntType>
	WaveletTreeRRR(const uIntType* src, size_type n, const Mapper& map) : n(n) {

	}

	/* member methods */
	/** test whether the i-th bit of val is set */
	bool test(unsigned int val, unsigned int i) const {
		assert (i < height);
		return (val & (1 << height - i - 1)) != 0;
	}

	/* utility methods */
private:
	/**
	 * build the basic fields from a given input Seq stroed in a BitStr
	 * @param bstr  encoded and bit-translated input
	 */
	void build_basic(const BitStrUint& bstr);

	/**
	 * recursively build Wavelet Tree bstrs of given level
	 */
	void build_level(unsigned int level, unsigned int length, uint offset);

	/* member fields */
private:
	unsigned int height = 0; /* height of Wavelet Tree */
	unsigned int max = 0; /* max value of input symbols */
	BitStrUint* bstrs; /* an array of BitStrings, with length height */
	size_t* OCC = nullptr; /* 0-based occurence of each symbol in alphabet, with length max + 2 */

	/* static methods */
public:
	template<typename uIntType>
	static uIntType max_value(const uIntType* src, size_t n) {
		uIntType maxV = 0;
		for(const uIntType* p = src; p != src + n; ++p)
			if(maxV < *p)
				maxV = *p;
		return maxV;
	}
};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* WAVELETTREERRR_H_ */
