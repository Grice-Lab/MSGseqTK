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
#include <algorithm>
#include <array>
#include "DNAalphabet.h"
#include "StringUtils.h"
#include "BAM.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::basic_string;
using std::istream;
using std::ostream;
using SAMtools::BAM;

/**
 * A DNASeq is a basic_string with nt16_t type values defined in DNAalphabet,
 * with added functions to handle unencoded strings directly
 * it stores the nt16_t values uncompressed, but support read/write to compressed manner
 */

class DNAseq: public std::basic_string<nt16_t> {
public:
	typedef std::array<size_t, DNAalphabet::SIZE> BASE_COUNT;
	/* nested types and enums */
	enum TRIM_MODE { FIVE_PRIME = 1, THREE_PRIME, BOTH_PRIME };

	/**
	 * Re-introduce base class methods
	 */
	using basic_string<nt16_t>::append;

	/* constructors */
	/** default constructor */
	DNAseq() = default;

	/** no virtual destructor since STL classes are not designed for inheritence */

	/** construct a DNAseq from copy of values */
	DNAseq(size_t n, value_type c) : basic_string(n, c)
	{   }

	/** constructing a DNAseq from a symbol string */
	explicit DNAseq(const string& str) {
		assign(str);
	}

	/** construting/converting from a basic_string<nt16_t> */
	DNAseq(const basic_string<nt16_t>& seq) : basic_string<nt16_t>(seq)
	{  }

	/** moving from a basic_string<nt16_t> */
	DNAseq(std::basic_string<nt16_t>&& seq) : basic_string<nt16_t>(seq)
	{  }

	/** destructor */
	virtual ~DNAseq() {  }

	/* member methods */
	/** copy assignment a DNAseq from a symbol string */
	DNAseq& operator=(const string& str);

	/** test whether position i is a valid */
	bool isValid(DNAseq::size_type i) const {
		return DNAalphabet::isValid((*this)[i]);
	}

	/** test whether the whole sequence is valid */
	bool isValid() const;

	/** test whether position i is a valid base (non-gap) */
	bool isBase(DNAseq::size_type i) const {
		return DNAalphabet::isBase((*this)[i]);
	}

	/** test whether the whole sequence is all base, no gap */
	bool isBase() const;

	/** test whether any base is gap */
	bool hasGap() const;

	/** reverse this object */
	DNAseq& reverse();

	/** Get a reverse copy */
	DNAseq reverse() const {
		DNAseq rSeq(*this);
		return rSeq.reverse();
	}

	/** Get a complement of this object */
	DNAseq& complement();

	/** Get a complement copy of */
	DNAseq complement() const {
		DNAseq cSeq(*this);
		return cSeq.complement();
	}

	/** reverse complement this object */
	DNAseq& revcom() {
		return reverse().complement();
	}

	/**
	 * Generate a reverse complement copy
	 */
	DNAseq revcom() const {
		DNAseq rcSeq(*this);
		return rcSeq.revcom();
	}

	/** get a substr copy of this object */
	DNAseq substr(size_t pos = 0, size_t len = npos) const {
		return basic_string<nt16_t>::substr(pos, len);
	}

	/**
	 * Get the decoded character at given position
	 * @param i  position
	 * @return  the decoded DNA character
	 */
	char decode(size_type i) const {
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
	using basic_string<nt16_t>::assign;

	/** Assign a string as a DNAseq */
	DNAseq& assign(const string& str);

	/**
	 * Append a seq string (symbols) to this DNAseq
	 * return the modified object
	 */
	DNAseq& append(const string& str);

	/** load a DNAseq from binary input */
	istream& load(istream& in);

	/** save a DNAseq to binary output */
	ostream& save(ostream& out) const;

	/** read a DNAseq from text input */
	istream& read(istream& in);

	/** write a DNAseq to text output */
	ostream& write(ostream& out) const {
		return out << toString();
	}

	/** encode this DNAseq to a BAM::seq_str in nt16 encoding */
	BAM::seq_str nt16Encode() const;

	/** load a DNAseq from binary input with compressed nt16 encoding */
	istream& nt16Load(istream& in);

	/** save a DNAseq to binary output with compressed nt16 encoding */
	ostream& nt16Save(ostream& out) const;

	/** remove invalid bases */
	DNAseq& removeInvalid();

	/** remove gap bases */
	DNAseq& removeGaps();

	/** get base count of a given region */
	BASE_COUNT baseCount(size_type start = 0, size_type len = npos) const;

	/* non-member functions */
	/** read a DNAseq from text input */
	friend istream& operator>>(istream& in, DNAseq& seq);

	/** write a DNAseq to text output */
	friend ostream& operator<<(ostream& out, const DNAseq& seq);

	/** concat two DNAseq */
	friend DNAseq operator+(const DNAseq& lhs, const DNAseq& rhs);

	/* static methods */
	/** decode a DNAseq from a nt16 bits encoded version and known length */
	static DNAseq nt16Decode(size_t L, const BAM::seq_str& seqNt16);

	/* static fields */
	static const DNAseq DNAgap; /* single base gap DNAseq */
};

inline ostream& operator<<(ostream& out, const DNAseq& seq) {
	return seq.write(out);
}

inline istream& operator>>(istream&in, DNAseq& seq) {
	return seq.read(in);
}

inline istream& DNAseq::load(istream& in) {
	return StringUtils::loadString(*this, in);
}

inline ostream& DNAseq::save(ostream& out) const {
	return StringUtils::saveString(*this, out);
}

inline bool DNAseq::isValid() const {
	return std::all_of(begin(), end(), [](DNAseq::value_type b) { return DNAalphabet::isValid(b); });
}

inline bool DNAseq::isBase() const {
	return std::all_of(begin(), end(), [](DNAseq::value_type b) { return DNAalphabet::isBase(b); });
}

inline bool DNAseq::hasGap() const {
	return std::any_of(begin(), end(), [](DNAseq::value_type b) { return DNAalphabet::isGap(b); });
}

inline DNAseq& DNAseq::operator=(const string& str) {
	assign(str);
	return *this;
}

inline DNAseq operator+(const DNAseq& lhs, const DNAseq& rhs) {
	DNAseq seq(lhs);
	seq += rhs;
	return seq;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNASEQ_H_ */
