/*
 * PrimarySeq.h
 *
 *  Created on: Jan 26, 2018
 *      Author: zhengqi
 */

#ifndef SRC_PRIMARYSEQ_H_
#define SRC_PRIMARYSEQ_H_

#include <string>
#include <iostream>
#include <cstdint>
#include "DNAalphabet.h"

namespace EGriceLab {
namespace MSGseqClean {

using std::basic_string;
using std::string;
using std::istream;
using std::ostream;

/**
 * A PrimarySeq contains a non-encoded DNA string with additional biological related information
 * such as name, description and quality
 */

class PrimarySeq {
	typedef basic_string<uint8_t> qstring;

public:
	/* constructors */
	/** default constructor */
	PrimarySeq() = default;

	/** constructing a PrimarySeq with given name and optional desc and encoded quality string */
	PrimarySeq(const string& seq, const string& name, const string& desc = "",
			const string& qStr = "", uint8_t qShift = DEFAULT_Q_SHIFT);

	/* getters and setters */
	const string& getSeq() const {
		return seq;
	}

	void setSeq(const string& seq) {
		this->seq = seq;
	}

	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	const string& getDesc() const {
		return desc;
	}

	void setDesc(const string& desc) {
		this->desc = desc;
	}

	qstring getQual() const {
		return !qual.empty() ? qual : qstring{ seq.length(), DEFAULT_Q_SCORE };
	}

	uint8_t getQShift() const {
		return qShift;
	}

	void setQShift(uint8_t qshift) {
		this->qShift = qShift;
	}

	void setQual(const qstring& qual) {
		this->qual = qual;
	}

	/**
	 * Append a seq string and quality scores to this PrimarySeq
	 * return the modified object
	 */
	PrimarySeq& append(const string& str, const string& qStr);

	/** load a DNAseq from binary input, override the old one */
	istream& load(istream& in);

	/** save a DNAseq to binary output */
	ostream& save(ostream& out) const;

	/** get qual get given position */
	uint8_t getQual(size_t i) const {
		return qual[i];
	}

	/** set qual at given position */
	void setQual(size_t i, uint8_t q) {
		qual[i] = q;
	}

	/** get encoded quality string */
	string getQStr() const;

	/** set qual using encoded string */
	void setQStr(const string& qStr);

	/** get a reverse-complement copy of this PrimarySeq object */
	PrimarySeq revcom() const;

	/* non-member functions */
	/** check equivenent of two PrimarySeq based on all fields */
	friend bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs);


	/* member fields */
private:
	string seq;
	string name;
	string desc;
	qstring qual;
	uint8_t qShift = DEFAULT_Q_SHIFT; // C++11

	/* static members */
public:
	static const uint8_t DEFAULT_Q_SCORE = 30;
	static const uint8_t DEFAULT_Q_SHIFT = 33;
};

inline bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return lhs.seq == rhs.seq && lhs.name == rhs.name && lhs.desc == rhs.desc && lhs.qual == rhs.qual;
}

inline bool operator!=(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return !(lhs == rhs);
}

inline PrimarySeq PrimarySeq::revcom() const {
	return PrimarySeq(DNAalphabet::revcom(seq), name, desc,
			string(qual.rbegin(), qual.rend()), qShift);
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_PRIMARYSEQ_H_ */
