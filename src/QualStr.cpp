/*
 * QualityString.cpp
 *
 *  Created on: Jun 14, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "QualStr.h"

namespace EGriceLab {
namespace MSGseqTK {
namespace quality {

QualStr encode(const string& qStr, uint8_t qShift) {
	QualStr qual(qStr.length(), 0);
	std::transform(qStr.begin(), qStr.end(), qual.begin(),
			[=](string::value_type q) { return q - qShift; });
	return qual;
}

string decode(const QualStr& qual, uint8_t qShift) {
	string qStr(qual.length(), '\0');
	std::transform(qual.begin(), qual.end(), qStr.begin(),
			[=](QualStr::value_type q) { return q + qShift; });
	return qStr;
}

istream& read(QualStr& qual, istream& in, uint8_t qShift) {
	string qStr;
	in >> qStr;
	qual = encode(qStr, qShift);
	return in;
}

} /* namespace quality */
} /* namespace MSGseqTK */
} /* namespace EGriceLab */
