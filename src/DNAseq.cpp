/*
 * DNAseq.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */
#include <algorithm>
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGSeqClean {

string DNAseq::decode() const {
	string seq;
	seq.reserve(length());
	for(uint8_t b : *this)
		seq.push_back(DNAalphabet::decode(b));
	return seq;
}

DNAseq DNAseq::revcom() const {
	DNAseq rcSeq;
	rcSeq.reserve(length());
	for(DNAseq::const_iterator b = rbegin(); b != rend(); ++b)
		rcSeq.push_back(DNAalphabet::complement(*b));
	return rcSeq;
}

DNAseq& DNAseq::append(const string& str) {
	for(char s : str)
		push_back(DNAalphabet::encode(s));
	return *this;
}

istream& operator>>(istream& in, DNAseq& seq) {
	string str;
	in >> str;
	seq.clear();
	seq.append(str);
	return in;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */


