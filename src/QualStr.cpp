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

QualStr& QualStr::assign(const string& str) {
	clear();
	for(char q : str)
		push_back(q - qShift);
	return *this;
}

string QualStr::decode() const {
	string qStr;
	qStr.reserve(length());
	for(QualStr::value_type q : *this)
		qStr.push_back(q + qShift);
	return qStr;
}

istream& operator>>(istream& in, QualStr& qual) {
	string str;
	in >> str;
	qual.assign(str);
	return in;
}

QualStr& QualStr::reverse() {
	std::reverse(begin(), end());
	return *this;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
