/*
 * DNAseq.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */
#include <algorithm>
#include <cassert>
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
	for(DNAseq::const_reverse_iterator b = rbegin(); b != rend(); ++b)
		rcSeq.push_back(DNAalphabet::complement(*b));
	return rcSeq;
}

DNAseq& DNAseq::assign(const string& str) {
	clear();
	reserve(str.length());
	for(char s : str)
		push_back(DNAalphabet::encode(s));
	return *this;
}

DNAseq& DNAseq::append(const string& str) {
	for(char s : str)
		push_back(DNAalphabet::encode(s));
	return *this;
}

istream& operator>>(istream& in, DNAseq& seq) {
	string str;
	in >> str;
	seq.assign(str);
	return in;
}

DNAseq& DNAseq::removeInvalid() {
	erase(std::remove(begin(), end(), 0), end());
	return *this;
}

DNAseq& DNAseq::removeGaps() {
	erase(std::remove(begin(), end(), DNAalphabet::N), end());
	return *this;
}

DNAseq& DNAseq::compressGaps(int minNGap) {
	assert(minNGap > 1);
	DNAseq::iterator gap_start, gap_end;
	while((gap_start = std::search_n(begin(), end(), minNGap, DNAalphabet::N)) != end()) { /* a consecutive gap found */
		for(gap_end = gap_start; *gap_end == DNAalphabet::N; ++gap_end) /* find gap end */
			continue;
		erase(gap_start + 1, gap_end);
	}
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */


