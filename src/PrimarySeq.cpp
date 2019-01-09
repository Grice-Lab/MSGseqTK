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
namespace MSGseqTK {

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
	qual.reverse();
	return *this;
}

PrimarySeq& PrimarySeq::complement() {
	seq.complement();
	return *this;
}

PrimarySeq& PrimarySeq::trunc(size_t start, size_t len) {
	seq = seq.substr(start, len);
	if(!qual.empty())
		qual = qual.substr(start, len);
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

