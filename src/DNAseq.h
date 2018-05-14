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
namespace MSGseqClean {

using std::string;
using std::istream;
using std::ostream;

/**
 * A DNASeq is a basic_string in int8_t type and always uses the DNAalphabet,
 * with added functions to handle unencoded strings directly
 */

class DNAseq: public std::basic_string<int8_t> {
public:
	/* constructors */
	/** default constructor */
	DNAseq() = default;

	/** virtual destructor */
	virtual ~DNAseq() {  }

	/** constructing a DNAseq from a symbol string */
	explicit DNAseq(const string& str) {
		assign(str);
	}

	/** additional copy assignment operator */
	DNAseq& operator=(const string& str);

	/* member methods */
	/**
	 * test whether position i is a valid
	 */
	bool isValid(DNAseq::size_type i) const {
		return DNAalphabet::isValid((*this)[i]);
	}

	/**
	 * test whether the whole sequence is valid
	 */
	bool allValid() const;

	/**
	 * test whether position i is a valid base (non-gap)
	 */
	bool isBase(DNAseq::size_type i) const {
		return DNAalphabet::isBase((*this)[i]);
	}

	/**
	 * test whether the whole sequence is all base, no gap
	 */
	bool allBase() const;


	/**
	 * Generate a reverse complement copy
	 * return a new copy in reverse complement version
	 */
	DNAseq revcom() const;

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

	/** Re-introduce all base class assign methods */
	using basic_string<int8_t>::assign;

	/** Assign a string as a DNAseq */
	DNAseq& assign(const string& str);

	/**
	 * Re-introduce all base class append methods
	 */
	using basic_string<int8_t>::append;

	/**
	 * Append a seq string (symbols) to this DNAseq
	 * return the modified object
	 */
	DNAseq& append(const string& str);

	/** load a DNAseq from binary input, override the old one */
	istream& load(istream& in);

	/** save a DNAseq to binary output */
	ostream& save(ostream& out) const;

	/** remove invalid bases (zero value) */
	DNAseq& removeInvalid();

	/** remove gap bases */
	DNAseq& removeGaps();

	/** remove genomic gaps (runs of Ns) and replace with a single N to save storage space */
	DNAseq& compressGaps(int minNGap = MIN_GENOME_GAP);

	/* non-member functions */
	/** read a DNAseq from text input, override the old one */
	friend istream& operator>>(istream& in, DNAseq& seq);

	/** write a DNAseq to text output */
	friend ostream& operator<<(ostream& out, const DNAseq& seq);

	/* static member fields */
	static const int MIN_GENOME_GAP = 10;

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

inline bool DNAseq::allValid() const {
	return std::all_of(begin(), end(), [](int8_t b) { return DNAalphabet::isValid(b); } );
}

inline bool DNAseq::allBase() const {
	return std::all_of(begin(), end(), [](int8_t b) { return DNAalphabet::isBase(b); } );
}

inline DNAseq& DNAseq::operator=(const string& str) {
	assign(str);
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNASEQ_H_ */
