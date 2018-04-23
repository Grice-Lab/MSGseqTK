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
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGSeqClean {

using std::string;
using std::istream;
using std::ostream;

/**
 * A PrimarySeq contains a DNAseq with additional biological related information,
 * such as name, description and quality, if available
 * PrimarySeq is not a DNAseq since it doesn't support many methods and operators such as append
 */

class PrimarySeq {
	typedef basic_string<uint8_t> qstring;

public:
	/* constructors */
	/** default constructor */
	PrimarySeq() = default;

	/** virtual destructor */
	virtual ~PrimarySeq() {  };

	/** constructing a PrimarySeq with given name and optional desc and encoded quality string */
	PrimarySeq(const string& seqStr, const string& name, const string& desc = "", const string& qStr = "", uint8_t qShift = DEFAULT_Q_SHIFT);

	/* getters and setters */
	const DNAseq& getSeq() const {
		return seq;
	}

	void setSeq(const DNAseq& seq) {
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

	const qstring& getQual() const {
		return qual;
	}

	void setQual(const qstring& qual) {
		this->qual = qual;
	}

	/**
	 * Re-introduce all base class append methods
	 */
	using DNAseq::append;

	/**
	 * Append a seq string and quality scores to this PrimarySeq
	 * return the modified object
	 */
	PrimarySeq& append(const string& str, const string& qStr);

	/** load a DNAseq from binary input, override the old one */
	istream& load(istream& in);

	/** save a DNAseq to binary output */
	ostream& save(ostream& out) const;

	/** get seq as string */
	string getSeqStr() const {
		return seq.toString();
	}

	/** set seq as given string, override old value */
	void setSeqStr(const string& str) {
		seq = str;
	}

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

	/* non-member functions */
	/** check equivenent of two PrimarySeq based on all fields */
	friend bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs);


	/* member fields */
private:
	DNAseq seq;
	string name;
	string desc;
	qstring qual;
	uint8_t qShift = DEFAULT_Q_SHIFT;

	/* static members */
public:
	static const uint8_t DEFAULT_Q_SCORE = 30;
	static const uint8_t DEFAULT_Q_SHIFT = 33;
};

inline bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return lhs.name == rhs.name && lhs.desc == rhs.desc && lhs.qual == rhs.qual &&
			dynamic_cast<const DNAseq&>(lhs) == dynamic_cast<const DNAseq&>(rhs);
}

inline bool operator!=(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_PRIMARYSEQ_H_ */
