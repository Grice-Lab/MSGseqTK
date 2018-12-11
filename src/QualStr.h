/*
 * QualStr.h
 *
 *  Created on: Jun 14, 2018
 *      Author: zhengqi
 */

#ifndef SRC_QUALSTR_H_
#define SRC_QUALSTR_H_

#include <string>
#include <cstdint>
#include <iostream>
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::istream;
using std::ostream;

/**
 * A Phred Quality String class
 */
class QualStr: public std::basic_string<uint8_t> {
public:
	/** default constructor */
	QualStr() = default;

	/** construct a quality string from a regular string */
	explicit QualStr(const string& str, uint8_t qShift = DEFAULT_Q_SHIFT)
	: qShift(qShift) {
		assign(str);
	}

	/** construct a QualStr with fixed number of values */
	explicit QualStr(size_type n, uint8_t val = DEFAULT_Q_SCORE)
	: std::basic_string<uint8_t>(n, val)
	 {   }

	/** additional copy assignment operator */
	QualStr& operator=(const string& str) {
		return assign(str);
	}

	/** no virtual destructor since STL classes are not designed for inheritence */

	/* member methods */
	/** get qShift */
	uint8_t getQShift() const {
		return qShift;
	}

	/** set qShift */
	void setQShift(uint8_t qShift) {
		this->qShift = qShift;
	}

	/** assign a string to this QualStr */
	QualStr& assign(const string& str);

	/** decode the quality as a std::string */
	string decode() const;

	/** get this QualStr as a string, alias to decode() */
	string toString() const {
		return decode();
	}

	/** test whether a given loc is valid */
	bool isValid(QualStr::size_type i) const {
		return (*this)[i] != INVALID_Q_SCORE;
	}

	/** test whether all quals are valid */
	bool isValid() const {
		return std::none_of(begin(), end(), [] (QualStr::value_type q) { return q == INVALID_Q_SCORE; });
	}

	/** get a substr copy of this QualStr */
	QualStr substr(size_t pos = 0, size_t len = npos) const;

	/** save this QualStr to binary output */
	ostream& save(ostream& out) const {
		return StringUtils::saveString(*this, out);
	}

	/** load a QualStr from binary input */
	istream& load(istream& in) {
		return StringUtils::loadString(*this, in);
	}

	/** write this QualStr to formatted output */
	ostream& write(ostream& out) const {
		return out << decode();
	}

	/** read a QualStr from formatted input */
	istream& read(istream& in);

	/** reverse this QualStr */
	QualStr& reverse();

	/** get a reversed copy of this QualStr */
	QualStr reverse() const {
		QualStr rQual(*this);
		return rQual.reverse();
	}

	/* non-member methods */
	/** formatted output */
	friend ostream& operator<<(ostream& out, const QualStr& qual);
	/** formatted input */
	friend istream& operator>>(istream& in, QualStr& qual);

	/* member fields */
private:
	uint8_t qShift = DEFAULT_Q_SHIFT;

public:
	/* static members */
	static const uint8_t DEFAULT_Q_SCORE = 30;
	static const uint8_t DEFAULT_Q_SHIFT = 33;
	static const uint8_t MIN_Q_SCORE = 2; /* prevent Inf */
	static const uint8_t INVALID_Q_SCORE = 0xFF;
};

inline ostream& operator<<(ostream& out, const QualStr& qual) {
	return qual.write(out);
}

inline istream& operator>>(istream& out, QualStr& qual) {
	return qual.read(out);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_QUALSTR_H_ */
