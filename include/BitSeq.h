/*
 * BitSeq.h
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#ifndef BITSEQ_H_
#define BITSEQ_H_
#include <cstddef>
#include <iostream>

namespace EGriceLab {
namespace libSDS {

using std::istream;
using std::ostream;

/**
 * A Succinct data structure for storing binary sequences with rank and select support
 * this is an abstracted class with many unimplemented methods
 */
class BitSeq {
public:
	/** destructor */
	virtual ~BitSeq() {  }

	/* member methods */
	/**
	 * Abstract method
	 * get the length of this BitSeq in bits
	 */
	virtual size_t length() const = 0;

	/**
	 * Abstract method
	 * get number of one bits (ones)
	 */
	virtual size_t ones() const = 0;

	/**
	 * Abstract method
	 * get the size of the structure in bytes
	 */
	virtual size_t getBytes() const = 0;

	/**
	 * Abstract method
	 * get # of ones until position i
	 * @param i  position
	 */
	virtual size_t rank1(size_t i) const = 0;

	/**
	 * get # of zeros until position i
	 */
	virtual size_t rank0(size_t i) const;

	/**
	 * Abstract method
	 * get the position of the i-th one
	 * @param r  required order/rank of 1
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > ones
	 */
	virtual size_t select1(size_t r) const = 0;

	/**
	 * Abstract method
	 * get the position of the i-th zero
	 * @param r  the order/rank of the 0
	 * @return  position in this BitSeq, or -1 if r = 0, or len if r > zeros
	 */
	virtual size_t select0(size_t r) const = 0;

	/**
	 * get the i-th bit of this BitSeq by testing the rank1 difference
	 * @param i  pos
	 * @return  true if it is i-th bit is significant (one)
	 */
	virtual bool access(size_t i) const;

	/* test bit i, alias to access */
	bool test(size_t i) const {
		return access(i);
	}

	/**
	 * get the next position starting from i that is one
	 * @param start  starting pos
	 * @return  next one position
	 */
	virtual size_t selectNext1(size_t start) const;

	/**
	 * get the next position starting from i that is zero
	 * @param start  starting pos
	 * @return  next zero position
	 */
	virtual size_t selectNext0(size_t start) const;

	/**
	 * get the previous position starting from i that is one
	 * @param start  starting pos
	 * @return  prev one position
	 */
	virtual size_t selectPrev1(size_t start) const;

	/**
	 * get the previous position starting from i that is zero
	 * @param start  starting pos
	 * @return  prev zero position
	 */
	virtual size_t selectPrev0(size_t start) const;

	/**
	 * abstract method
	 * save this BitSeq to binary output
	 */
	virtual ostream& save(ostream& out) const = 0;

	/**
	 * abstract method
	 * load data from a binary input
	 */
	virtual istream& load(istream& in) = 0;
};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSEQ_H_ */
