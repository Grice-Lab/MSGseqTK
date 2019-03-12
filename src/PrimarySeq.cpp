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
	dna::load(seq, in);
	StringUtils::loadString(name, in);
	StringUtils::loadString(desc, in);
	quality::load(qual, in);
	return in;
}

ostream& PrimarySeq::save(ostream& out) const {
	dna::save(seq, out);
	StringUtils::saveString(name, out);
	StringUtils::saveString(desc, out);
	quality::save(qual, out);
	return out;
}

PrimarySeq& PrimarySeq::reverse() {
	dna::reverse(seq);
	quality::reverse(qual);
	return *this;
}

PrimarySeq& PrimarySeq::complement() {
	dna::complement(seq);
	return *this;
}

PrimarySeq& PrimarySeq::trunc(size_t start, size_t len) {
	seq = seq.substr(start, len);
	if(!qual.empty())
		qual = qual.substr(start, len);
	return *this;
}

PrimarySeq& PrimarySeq::trimNameExt(size_t len) {
	if(name.length() > len) // anything left-over
		name.erase(name.length() - len);
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

