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
	for(char q : str) {
		uint8_t qScore = q - qShift;
		if(qScore < MIN_Q_SCORE)
			qScore = MIN_Q_SCORE;
		push_back(qScore);
	}
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

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
