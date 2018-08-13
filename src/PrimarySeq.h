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
#include "MSGseqTKConst.h"
#include "QualStr.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::basic_string;
using std::string;
using std::istream;
using std::ostream;

/**
 * A PrimarySeq contains a non-encoded DNA string with additional biological related information
 * such as name, description and quality
 */

class PrimarySeq {
public:
	/* constructors */
	/** default constructor */
	PrimarySeq() = default;

	/** construct a PrimarySeq with given seq, qual, name and optional desc and encoded quality string */
	PrimarySeq(const string& seqStr, const string& name, const string& desc = "",
			const string& qStr = "", uint8_t qShift = QualStr::DEFAULT_Q_SHIFT);

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

	QualStr getQual() const;

	void setQual(const QualStr& qual) {
		this->qual = qual;
	}

	size_t length() const {
		return seq.length();
	}

	bool empty() const {
		return seq.empty();
	}

	/** load a DNAseq from binary input, override the old one */
	istream& load(istream& in);

	/** save a DNAseq to binary output */
	ostream& save(ostream& out) const;

	/** reverse this Primary seq */
	PrimarySeq& reverse();

	/** get a reverse copy of this Primary seq */
	PrimarySeq reverse() const {
		PrimarySeq rSeq(*this);
		return rSeq.reverse();
	}

	/** complement this Primary Seq */
	PrimarySeq& complement();

	/** get a complement copy of this Primary Seq */
	PrimarySeq complement() const {
		PrimarySeq cSeq(*this);
		return cSeq.complement();
	}

	/** reverse-complement this object */
	PrimarySeq& revcom() {
		return reverse().complement();
	}

	/** get a reverse-complement copy of this PrimarySeq object */
	PrimarySeq revcom() const {
		PrimarySeq rcSeq(*this);
		return rcSeq.revcom();
	}

	/* non-member functions */
	/** check equivenent of two PrimarySeq based on all fields */
	friend bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs);

	/* member fields */
private:
	DNAseq seq;
	string name;
	string desc;
	QualStr qual;
};

inline bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return lhs.seq == rhs.seq && lhs.name == rhs.name && lhs.desc == rhs.desc && lhs.qual == rhs.qual;
}

inline bool operator!=(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_PRIMARYSEQ_H_ */
