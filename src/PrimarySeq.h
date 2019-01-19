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
#include <algorithm>
#include <cctype>
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

	/** construct a PrimarySeq with given data */
	PrimarySeq(const DNAseq& seq, const string& name, const string& desc,
			const QualStr& qual, uint8_t qShift = QualStr::DEFAULT_Q_SHIFT)
	: seq(seq), name(name), desc(desc), qual(qual)
	{
		this->qual.setQShift(qShift);
		fixQual();
	}

	/** delegating construct a PrimarySeq with given seq, qual, name and optional desc and encoded quality string */
	PrimarySeq(const DNAseq& seq, const string& name, const string& desc = "",
			const string& qStr = "", uint8_t qShift = QualStr::DEFAULT_Q_SHIFT)
	: PrimarySeq(seq, name, desc, QualStr(qStr, qShift), qShift)
	{  }

	/** delegating construct using a seq string instead of DNAseq */
	PrimarySeq(const string& seqStr, const string& name, const string& desc = "",
			const string& qStr = "", uint8_t qShift = QualStr::DEFAULT_Q_SHIFT)
	: PrimarySeq(DNAseq(seqStr), name, desc, QualStr(qStr, qShift), qShift)
	{  }

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

	bool hasQual() const {
		return qual.length() == seq.length();
	}

	const QualStr& getQual() const {
		return qual;
	}

	void setQual(const QualStr& qual) {
		this->qual = qual;
	}

	size_t length() const {
		return seq.length();
	}

	bool empty() const {
		return seq.empty();
	}

	/** check and fix QualStr */
	PrimarySeq& fixQual() {
		if(qual.length() != seq.length())
			qual.resize(seq.length(), QualStr::DEFAULT_Q_SCORE);
		return *this;
	}

	DNAseq::value_type getBase(DNAseq::size_type i) const {
		return seq[i];
	}

	QualStr::value_type getQvalue(QualStr::size_type i) const {
		return qual[i];
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

	/**
	 * truncate this PrimarySeq to a substring
	 * @param start  0-based start of substring
	 * @param len  length of substring
	 */
	PrimarySeq& trunc(size_t start, size_t len);

	/**
	 * get a substring copy of this PrimarySeq
	 */
	PrimarySeq substr(size_t start, size_t len) const {
		return qual.empty() ?
				PrimarySeq(seq.substr(start, len), name, desc, qual) :
				PrimarySeq(seq.substr(start, len), name, desc, qual.substr(start, len));
	}

	/** trim name extension */
	PrimarySeq& trimNameExt(size_t len = ILLUMINA_PAIR_ID_EXT_LEN);

	/* non-member functions */
	/** check equivenent of two PrimarySeq based on all fields */
	friend bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs);

	/* member fields */
private:
	DNAseq seq;
	string name;
	string desc;
	QualStr qual;

	/* static fields */
public:
	static const char FILL_CHAR = '.';
	static const size_t ILLUMINA_PAIR_ID_EXT_LEN = 2; // potential illumina paired-read ID ext length (/1/2 or .1.2)

	/* static methods */
public:
	/** format a given seq name for display to remove white space characters */
	static string formatName(string name) {
		std::replace_if(name.begin(), name.end(),
				[] (char c) { return ::isspace(c); },
				FILL_CHAR);
		return name;
	}
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
