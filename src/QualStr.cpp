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

const double QualStr::PHRED_SCALE = -10;

QualStr& QualStr::assign(const string& str) {
	resize(str.length());
	std::transform(str.begin(), str.end(), begin(),
			[=](string::value_type q) { return std::max<uint8_t>(q - qShift, MIN_Q_SCORE); });
	return *this;
}

string QualStr::decode() const {
	string qStr(length(), '\0');
	std::transform(begin(), end(), qStr.begin(),
			[=](value_type q) { return q + qShift; });
	return qStr;
}

QualStr QualStr::substr(size_t pos, size_t len) const {
	QualStr seg;
	if(pos + len >= length())
		len = length() - pos;
	seg.resize(len);
	std::copy_n(begin() + pos, len, seg.begin());
	return seg;
}

istream& QualStr::read(istream& in) {
	string str;
	in >> str;
	assign(str);
	return in;
}

QualStr& QualStr::reverse() {
	std::reverse(begin(), end());
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
