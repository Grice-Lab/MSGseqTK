/*
 * DNASeq.h
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */

#ifndef SRC_DNASEQ_H_
#define SRC_DNASEQ_H_

#include <string>
#include <iostream>
#include "DNAalphabet.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGSeqClean {

using std::string;
using std::istream;
using std::ostream;

/**
 * A DNASeq is a basic_string in uint8_t type and always uses the DNAalphabet,
 * with added functions to handle unencoded strings directly
 */

class DNAseq: public std::basic_string<uint8_t> {
public:
	/* constructors */
	/** default constructor */
	DNAseq() = default;

	/** virtual destructor */
	virtual ~DNAseq() {  }

	/** constructing a DNAseq from a symbol string */
	explicit DNAseq(const string& str) {
		append(str);
	}

	/** default copy assignment constructor */
	DNAseq& operator=(const DNAseq& other) = default;

	/** copy assignment constructor from a encoded string */
	DNAseq& operator=(const string& str) {
		clear();
		append(str);
		return *this;
	}

	/* member methods */
	/**
	 * test whether position i is a valid
	 */
	bool isValid(DNAseq::size_type i) const {
		return (*this)[i] > 0;
	}

	/**
	 * Generate a reverse complement copy
	 * return a new copy in reverse complement version
	 */
	DNAseq revcom() const noexcept;

	/**
	 * Get the decoded character at given position
	 * @param i  position
	 * @return  the decoded DNA character
	 */
	char decode(DNAseq::size_type i) const {
		return DNAalphabet::decode((*this)[i]);
	}

	/**
	 * decode the entire seq to symbols
	 */
	string decode() const;

	/**
	 * Alias of decode()
	 */
	string toString() const {
		return decode();
	}

	/**
	 * Re-introduce all base class append methods
	 */
	using basic_string<uint8_t>::append;

	/**
	 * Append a seq string (symbols) to this DNAseq
	 * return the modified object
	 */
	DNAseq& append(const string& str);

	/** load a DNAseq from binary input, override the old one */
	istream& load(istream& in);

	/** save a DNAseq to binary output */
	ostream& save(ostream& out) const;

	/* non-member functions */
	/** read a DNAseq from text input, override the old one */
	friend istream& operator>>(istream& in, DNAseq& seq);

	/** write a DNAseq to text output */
	friend ostream& operator<<(ostream& out, const DNAseq& seq);

};

inline ostream& operator<<(ostream& out, const DNAseq& seq) {
	return out << seq.decode();
}

inline istream& DNAseq::load(istream& in) {
	clear();
	return StringUtils::loadString(*this, in);
}

inline ostream& DNAseq::save(ostream& out) const {
	return StringUtils::saveString(*this, out);
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNASEQ_H_ */
