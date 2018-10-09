/*
 * BitSeqGGMN.h
 *
 *  Created on: Sep 24, 2018
 *      Author: zhengqi
 */

#ifndef BITSEQGGMN_H_
#define BITSEQGGMN_H_

#include <algorithm>
#include "BitSeq.h"
#include "BitStr.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

/**
 * An uncompressed implementation of BitSeq using the algorithm described in
 * [1] R. González, Sz. Grabowski, V. Mäkinen, and G. Navarro. Practical implementation of
 * rank and select queries. In Posters WEA, pages 27–38, 2005
 * and also described in
 * [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 */
class BitSeqGGMN: public BitSeq {
public:
	typedef uint data_type;
	/* constructors */
	/** default constructor */
	BitSeqGGMN() = default;

	/** copy constructor */
	BitSeqGGMN(const BitSeqGGMN& other);

	/**
	 * constructing from a BitStr of suitable type and factor rate
	 * @param bstr  BitStr of uint
	 * @param super-block factor, if 0, automatic determined
	 */
	BitSeqGGMN(const BitStr<uint>& bstr, size_t factor = 0);

	/**
	 * constructing from a BitStr of suitable type and factor rate
	 * @param bstr  BitStr of type T
	 * @param super-block factor, if 0, automatic determined
	 */
	template<typename oIntType>
	BitSeqGGMN(const BitStr<oIntType>& bstr, size_t factor = 0) : bstr(bstr) {
		if(factor == 0)
			factor = bits(length() - 1);
		this->factor = factor;
		b = sizeof(data_type) * Wb;
		s = b * factor;
		buildRank();
	}

	/** virtual destructor */
	virtual ~BitSeqGGMN()
	{  }

	/* member methods */
	/**
	 * get the length of this BitSeq in bits
	 * @override  base class virtual method
	 */
	virtual size_t length() const {
		return bstr.length();
	}

	/**
	 * get number of one bits (ones)
	 * @override  base class virtual method
	 */
	virtual size_t ones() const {
		return bstr.count();
	}

	/**
	 * get the size of the structure in bytes
	 * @override  base class virtual method
	 */
	virtual size_t getBytes() const;

	/**
	 * Abstract method
	 * get # of ones until position i
	 * @param i  position
	 * @override  base class virtual method
	 * @return  number of ones until position i,
	 * use level-one binary-search, level-2 and level-3 popcount and bit search
	 */
	virtual size_t rank1(size_t i) const;

	/**
	 * get the position of the i-th one
	 * @param r  the order/rank of the 1
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > ones
	 * @override  base class virtual method
	 * first binary search over first level rank structure
	 * then sequential search using popcount over a int
	 * then sequential search using popcount over a char
	 * then sequential search bit by bit
	 */
	virtual size_t select1(size_t r) const;

	/**
	 * get the position of the i-th zero
	 * @param r  the order/rank of the 0
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > zeros
	 * @override  base class virtual method
	 * first binary search over first level rank structure
	 * then sequential search using popcount over a int
	 * then sequential search using popcount over a char
	 * then sequential search bit by bit
	 */
	virtual size_t select0(size_t r) const;

	/**
	 * get the i-th bit of this BitSeq by testing the rank1 difference
	 * @param i  pos
	 * @return  true if it is i-th bit is significant (one)
	 * @override  base class virtual method
	 */
	virtual bool access(size_t i) const {
		return bstr.get(i);
	}

	/**
	 * get the next position starting from i that is one
	 * @param start  starting pos
	 * @return  next one position
	 * @override  base class virtual method
	 * this method is optimized for uint 32bits machines
	 */
	virtual size_t selectNext1(size_t start) const;

	/**
	 * get the previous position starting from i that is one
	 * @param start  starting pos
	 * @return  prev one position
	 * @override  base class virtual method
	 * this method is optimized for uint 32bits machines
	 */
	virtual size_t selectPrev1(size_t start) const;

	/**
	 * save this BitSeq to binary output
	 * @param out  binary output
	 * @override  base class virtual method
	 */
	virtual ostream& save(ostream& out) const;

	/**
	 * load data from a binary input
	 * @param in  binary input
	 * @override  base class virtual method
	 */
	virtual istream& load(istream& in);

	/** get number of superblocks, deduced from numBits() */
	size_t numSuperBlocks() const {
		return (bstr.length() + s - 1) / s;
	}

	/** swap this BitSeqGGMN with another object */
	void swap(BitSeqGGMN& other) {
		std::swap(bstr, other.bstr);
		std::swap(Rs, other.Rs);
		std::swap(factor, other.factor);
		std::swap(b, other.b);
		std::swap(s, other.s);
	}

	/* internal help methods */
private:

	/**
	 * build rank index on super-block level
	 */
	void buildRank();

	/**
	 * Build rank on give block
	 * @param start  start of sub-block
	 * @param len  length of sub-block
	 * @return  number of on bits (ones) in this sub-block
	 */
	size_t buildRank(size_t start, size_t len);

	/* member fields */
private:
	BitStr<data_type> bstr; /* default data type is uint for efficiency */
	size_t* Rs = nullptr; /* 1-based super-block rank array, storing the on bits in each block */
	size_t factor = 0;
	size_t b = 0; /* block size in bits */
	size_t s = 0; /* super-block size in bits per block */
};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSEQGGMN_H_ */
