/*
 * PrimarySeq.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: zhengqi
 */
#include <algorithm>
#include "PrimarySeq.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqClean {

PrimarySeq::PrimarySeq(const string& seq, const string& name, const string& desc,
		const string& qStr, uint8_t qShift)
: seq(seq), name(name), desc(desc), qual(qStr, qShift)
{   }

QualStr PrimarySeq::getQual() const {
	if(qual.length() == seq.length())
		return qual;
	else {
		return QualStr(seq.length(), QualStr::DEFAULT_Q_SCORE);
	}
}

istream& PrimarySeq::load(istream& in) {
	seq.load(in);
	StringUtils::loadString(name, in);
	StringUtils::loadString(desc, in);
	qual.load(in);
	return in;
}

ostream& PrimarySeq::save(ostream& out) const {
	seq.save(out);
	StringUtils::saveString(name, out);
	StringUtils::saveString(desc, out);
	qual.save(out);
	return out;
}

PrimarySeq& PrimarySeq::reverse() {
	seq.reverse();
	std::reverse(qual.begin(), qual.end());
	return *this;
}

PrimarySeq& PrimarySeq::complement() {
	seq.complement();
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

