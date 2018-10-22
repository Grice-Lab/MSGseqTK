/*
 * BitSeqGGMN.h
 *
 *  Created on: Sep 24, 2018
 *      Author: zhengqi
 */

#ifndef BITSEQGGMN_H_
#define BITSEQGGMN_H_

#include <vector>
#include <algorithm>
#include "BitSeq.h"
#include "BitStr.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

using std::vector;

/**
 * An uncompressed implementation of BitSeq using the algorithm described in
 * [1] R. González, Sz. Grabowski, V. Mäkinen, and G. Navarro. Practical implementation of
 * rank and select queries. In Posters WEA, pages 27–38, 2005
 * and also described in
 * [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 */
class BitSeqGGMN: public BitSeq {
public:
	typedef BitStr32::value_type data_type;
	/* constructors */
	/** default constructor */
	BitSeqGGMN() = default;

	/**
	 * constructing from a BitStr of suitable type and factor rate
	 * @param bstr  BitStr of uint
	 * @param super-block factor, if 0, automatic determined
	 */
	BitSeqGGMN(const BitStr32& bstr, size_t factor = 0);

	/**
	 * constructing from a BitStr of suitable type and factor rate
	 * @param bstr  BitStr of type T
	 * @param super-block factor, if 0, automatic determined
	 */
	template<typename oIntType>
	BitSeqGGMN(const BitStr<oIntType>& bstr, size_t factor = 0) : bstr(bstr) {
		n = bstr.length();
		ones = bstr.count();
		if(factor == 0)
			factor = bits(n - 1);
		this->factor = factor;
		b = sizeof(data_type) * Wb;
		s = b * factor;
		nRs = numSuperBlocks() + 1;
		buildRank();
	}

	/** virtual destructor */
	virtual ~BitSeqGGMN()
	{  }

	/**
	 * get the size of the structure in bytes
	 * @override  base class virtual method
	 */
	virtual size_t getBytes() const;

	/**
	 * get # of ones until position i (0-based, inclusive)
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

	/* import base class access methods */
	using BitSeq::access;

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
	 * @return  next one position, or n if not found
	 * @override  base class virtual method
	 * this method is optimized for 32 bits machines
	 */
	virtual size_t selectNext1(size_t start) const;

	/**
	 * get the previous position starting from i that is one
	 * @param start  starting pos
	 * @return  prev one position, or 0 if not found
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
		return (n + s - 1) / s;
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

	/* non-member methods */
	/* relationship operators */
public:
	friend bool operator==(const BitSeqGGMN& lhs, const BitSeqGGMN& rhs);

	/* member fields */
private:
	BitStr32 bstr; /* use BitStr32 for efficiency */
	size_t nRs = 0; /* size of Rs */
	vector<size_t> Rs; /* 1-based super-block rank vector, storing the on bits in each block, with Rs[0] = 0 */
	size_t factor = 0;
	size_t b = 0; /* block size in bits */
	size_t s = 0; /* super-block size in bits per block */
};

inline bool operator!=(const BitSeqGGMN& lhs, const BitSeqGGMN& rhs) {
	return !(lhs == rhs);
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSEQGGMN_H_ */
