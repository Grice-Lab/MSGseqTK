/*
 * PrimarySeq.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: zhengqi
 */

#include "PrimarySeq.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGSeqClean {

PrimarySeq::PrimarySeq(const string& seqStr, const string& name, const string& desc, const string& qStr, uint8_t qShift)
: seq(seqStr), name(name), desc(desc), qShift(qShift)
{
	for(char qCh : qStr)
		qual.push_back(qCh - qShift);
}

istream& PrimarySeq::load(istream& in) {
	/* load seq */
	seq.load(in);
	/* load other info */
	StringUtils::loadString(name, in);
	StringUtils::loadString(desc, in);
	StringUtils::loadString(qual, in);
	in.read((char *) &qShift, sizeof(uint8_t));

	return in;
}

ostream& PrimarySeq::save(ostream& out) const {
	/* save seq */
	seq.save(out);
	/* save other info */
	StringUtils::saveString(name, out);
	StringUtils::saveString(desc, out);
	StringUtils::saveString(qual, out);
	out.write((const char*) &qShift, sizeof(uint8_t));

	return out;
}

string PrimarySeq::getQStr() const {
	string qStr;
	qStr.reserve(length());
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

