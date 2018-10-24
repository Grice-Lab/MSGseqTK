/*
 * Seq.h
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#ifndef SEQ_H_
#define SEQ_H_
#include <cstddef>
#include <iostream>

namespace EGriceLab {
namespace libSDS {

using std::istream;
using std::ostream;

/**
 * A Succinct data structure for storing sequences in arbitary alphabet with rank and select support
 * this is an abstracted class with unimplemented methods
 */
class Seq {
public:
	typedef size_t size_type;
	/* constructors */
	/** default constructor */
	Seq() = default;

	/** destructor */
	virtual ~Seq()
	{ 	}

	/** basic constructor from given values */
	Seq(size_t n, size_t sigma) : n(n), sigma(sigma)
	{  }

	/* member methods */
	/**
	 * get the length of this BitSeq in bits
	 */
	virtual size_t length() const {
		return n;
	}

	/**
	 * get size of alphabet
	 */
	virtual size_t alphabetSize() const {
		return sigma;
	}

	/**
	 * get the size of the structure in bytes
	 */
	virtual size_t getBytes() const {
		return sizeof(n) + sizeof(sigma) + sizeof(this);
	}

	/**
	 * Abstract method
	 * get the i-th symbol in this Seq
	 */
	virtual size_t access(size_t i) const = 0;

	/**
	 * Retrieve the symbol at position i and its rank.
	 */
	virtual size_t access(size_t i, size_t & r) const;

	/**
	 * get #occurrence of symbol s until position i (0-based, inclusive)
	 * this naive implementation is based on counting using the abstract access
	 * @param s  symbol
	 * @param i  position
	 * @return rank of s till position i
	 */
	virtual size_t rank(size_t s, size_t i) const;

	/**
	 * get the position of the r-th occurrence of symbol s
	 * this naive impelmentation go through all symbols until found it
	 * @param s  symbol to search
	 * @param r  required occurrence/rank of s
	 * @return  position in this Seq, or -1 if i = 0, or len if i > total # of s
	 */
	virtual size_t select(size_t s, size_t r) const;

	/**
	 * get the next position of the r-th occurrence of symbol s
	 * this basic impelmentation depends on select and rank
	 * @param s  symbol to search
	 * @param start  start positin to search
	 * @return  next position in this Seq, or 0 if not found
	 */
	virtual size_t selectNext(size_t s, size_t start) const;

	/**
	 * get the max value of symbols in this Seq
	 */
	virtual size_t maxValue() const;

	/**
	 * save this BitSeq to binary output
	 */
	virtual ostream& save(ostream& out) const;

	/**
	 * load data from a binary input
	 */
	virtual istream& load(istream& in);

	/* non-member methods */
	/* relationship operators */
	friend bool operator==(const Seq& lhs, const Seq& rhs);

	/* member fields */
protected:
	size_t n = 0;     /* length of Seq */
	size_t sigma = 0; /* size of alphabet */
};

inline bool operator==(const Seq& lhs, const Seq& rhs) {
	return lhs.n == rhs.n && lhs.sigma == rhs.sigma;
}

inline bool operator!=(const Seq& lhs, const Seq& rhs) {
	return !(lhs == rhs);
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* SEQ_H_ */
