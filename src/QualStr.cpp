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
	clear();
	for(char q : str)
		push_back(std::max<uint8_t>(q - qShift, MIN_Q_SCORE));
	return *this;
}

string QualStr::decode() const {
	string qStr;
	qStr.reserve(length());
	for(QualStr::value_type q : *this)
		qStr.push_back(q + qShift);
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
