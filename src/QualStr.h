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
#include <cmath>
#include <algorithm>
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::istream;
using std::ostream;

/**
 * A Phred Quality String is just a convenient typedef with stand-alone functions
 */
typedef std::basic_string<uint8_t> QualStr;

namespace quality {
const uint8_t DEFAULT_Q_SCORE = 30;
const uint8_t DEFAULT_Q_SHIFT = 33;
const uint8_t MIN_Q_SCORE = 2; /* prevent Inf */
const uint8_t INVALID_Q_SCORE = 0xFF;
const uint8_t MAX_Q_SCORE = 250;
const double PHRED_SCALE = -10;

/** encode a QualStr from a string */
QualStr encode(const string& qStr, uint8_t qShift = DEFAULT_Q_SHIFT);

/** decode a QualStr to a string */
string decode(const QualStr& qual, uint8_t qShift = DEFAULT_Q_SHIFT);

/** alias to decode() */
inline
string toString(const QualStr& qual, uint8_t qShift = DEFAULT_Q_SHIFT) {
	return decode(qual, qShift);
}

/** save a QualStr to binary output */
inline
ostream& save(const QualStr& qual, ostream& out) {
	return StringUtils::saveString(qual, out);
}

/** load a QualStr from binary input */
inline
istream& load(QualStr& qual, istream& in) {
	return StringUtils::loadString(qual, in);
}

/** write this QualStr to formatted output */
inline
ostream& write(const QualStr& qual, ostream& out, uint8_t qShift = DEFAULT_Q_SHIFT) {
	return out << decode(qual, qShift);
}

/** read a QualStr from formatted input */
istream& read(QualStr& qual, istream& in, uint8_t qShift = DEFAULT_Q_SHIFT);

/** reverse a QualStr */
inline
QualStr& reverse(QualStr& qual) {
	std::reverse(qual.begin(), qual.end());
	return qual;
}

/** get a reversed copy of a QualStr */
inline
QualStr reverse(const QualStr& qual) {
	QualStr rQual(qual);
	return reverse(rQual);
}

/** formatted output */
inline
ostream& operator<<(ostream& out, const QualStr& qual) {
	return write(qual, out);
}

/** formatted input */
inline
istream& operator>>(istream& in, QualStr& qual) {
	return read(qual, in);
}

/* quality score transforming functions */
inline
double phredQ2P(double q) {
	return ::pow(10.0, q / PHRED_SCALE);
}

inline
double phredP2Q(double p) {
	return PHRED_SCALE * ::log10(p);
}

} /* namespace quality */
} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_QUALSTR_H_ */
