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
	explicit QualStr(size_type n = 0, uint8_t val = DEFAULT_Q_SHIFT)
	: std::basic_string<uint8_t>(n, val)
	 {   }

	/** additional copy assignment operator */
	QualStr& operator=(const string& str) {
		return assign(str);
	}

	virtual ~QualStr() {
		// TODO Auto-generated destructor stub
	}

	/* member methods */
	/** assign a string to this QualStr */
	QualStr& assign(const string& str);

	/** decode the quality as a std::string */
	string decode() const;

	/** get this QualStr as a string, alias to decode() */
	string toString() const {
		return decode();
	}

	/** save this QualStr to binary output */
	ostream& save(ostream& out) const {
		return StringUtils::saveString(*this, out);
	}

	/** load a QualStr from binary input */
	istream& load(istream& in) {
		return StringUtils::loadString(*this, in);
	}

	/** reverse this QualStr */
	QualStr& reverse();

	/** get a reversed copy of this QualStr */
	QualStr reverse() const {
		QualStr rQual(*this);
		return rQual.reverse();
	}

	/* non-member methods */
	friend ostream& operator<<(ostream& out, const QualStr& qual);

	friend istream& operator>>(istream& in, QualStr& qual);

	/* member fields */
private:
	uint8_t qShift = DEFAULT_Q_SHIFT;

public:
	/* static members */
	static const uint8_t DEFAULT_Q_SCORE = 30;
	static const uint8_t DEFAULT_Q_SHIFT = 33;
};

inline ostream& operator<<(ostream& out, const QualStr& qual) {
	return out << qual.decode();
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_QUALSTR_H_ */
