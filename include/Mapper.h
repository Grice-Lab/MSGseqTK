/*
 * Mapper.h
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#ifndef MAPPER_H_
#define MAPPER_H_
#include <cstddef>
#include <iostream>

namespace EGriceLab {
namespace libSDS {

using std::istream;
using std::ostream;

/* a basic mapper to map arbitary size of alphabet to arbitary size of alphabet */
class Mapper {
public:
	/** destructor */
	virtual ~Mapper()
	{  }

	/* member methods */
	/**
	 * Abstract method
	 * map/encode a character in input alphabet to output alphabet
	 */
	virtual unsigned int encode(unsigned int s) const = 0;

	/**
	 * Abstract method
	 * unmap/decode a character from output alphabet to input alphabet
	 */
	virtual unsigned int decode(unsigned int s) const = 0;

	/**
	 * get the size of the structure in bytes
	 */
	virtual size_t getBytes() const {
		return sizeof(this);
	}

	/**
	 * Abstract method
	 * save this object to binary output
	 */
	virtual ostream& save(ostream& out) const = 0;

	/**
	 * Abstract method
	 * load object from a binary input
	 */
	virtual istream& load(istream& in) = 0;
};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* MAPPER_H_ */
