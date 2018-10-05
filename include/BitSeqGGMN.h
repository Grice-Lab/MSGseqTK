/*
 * BitSeqGGMN.h
 *
 *  Created on: Sep 24, 2018
 *      Author: zhengqi
 */

#ifndef BITSEQGGMN_H_
#define BITSEQGGMN_H_

#include "BitSeq.h"
#include "BitStr.h"

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
	/* constructors */
	/** default constructor */
	BitSeqGGMN() = default;

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
	BitSeqGGMN(const BitStr<class oIntType>& bstr, size_t factor = 0);

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
	virtual size_t getSize() const;

	/**
	 * Abstract method
	 * get # of ones until position i
	 * @param i  position
	 * @override  base class virtual method
	 */
	virtual size_t rank1(const size_t i) const;

	/**
	 * get # of zeros until position i
	 * @param i  position
	 * @override  base class virtual method
	 */
	virtual size_t rank0(const size_t i) const;

	/**
	 * get the position of the i-th one
	 * @param i  the order of the 1
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > ones
	 * @override  base class virtual method
	 */
	virtual size_t select1(const size_t i) const;

	/**
	 * get the position of the i-th zero
	 * @param i  the order of the 0
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > zeros
	 * @override  base class virtual method
	 */
	virtual size_t select0(const size_t i) const;

	/**
	 * get the i-th bit of this BitSeq by testing the rank1 difference
	 * @param i  pos
	 * @return  true if it is i-th bit is significant (one)
	 * @override  base class virtual method
	 */
	virtual bool access(const size_t i) const;

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

	/* internal help methods */
private:
	/**
	 * build rank index on super-block level
	 */
	void BuildRank();

	/**
	 * Build rank on give block
	 * @param start  start of sub-block
	 * @param len  length of sub-block
	 * @return  number of on bits (ones) in this sub-block
	 */
	size_t BuildRankSub(size_t start, size_t len);

	/* member fields */
private:
	BitStr<uint> bstr; /* use uint size for efficiency */
	size_t factor = 0;
	size_t b = 0; /* block size per element */
	size_t s = 0; /* super-block size per element */
};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSEQGGMN_H_ */
