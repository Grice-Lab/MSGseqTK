/*
 * PrimarySeq.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: zhengqi
 */

#include "PrimarySeq.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqClean {

PrimarySeq::PrimarySeq(const string& seq, const string& name, const string& desc,
		const string& qStr, uint8_t qShift)
: seq(seq), name(name), desc(desc), qShift(qShift)
{
	for(char q : qStr)
		qual.push_back(q - qShift);
}

istream& PrimarySeq::load(istream& in) {
	StringUtils::loadString(seq, in);
	StringUtils::loadString(name, in);
	StringUtils::loadString(desc, in);
	StringUtils::loadString(qual, in);
	in.read((char *) &qShift, sizeof(uint8_t));

	return in;
}

ostream& PrimarySeq::save(ostream& out) const {
	StringUtils::saveString(seq, out);
	StringUtils::saveString(name, out);
	StringUtils::saveString(desc, out);
	StringUtils::saveString(qual, out);
	out.write((const char*) &qShift, sizeof(uint8_t));

	return out;
}

string PrimarySeq::getQStr() const {
	if(qual.empty())
		return string(seq.length(), DEFAULT_Q_SCORE + qShift);

	string qStr;
	qStr.reserve(seq.length());
	for(uint8_t q : qual)
		qStr.push_back(q + qShift);
	return qStr;
}

void PrimarySeq::setQStr(const string& qStr) {
	qual.clear();
	qual.reserve(qStr.length());
	for(char qCh : qStr)
		qual.push_back(qCh - qShift);
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

