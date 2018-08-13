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
namespace MSGseqTK {

string DNAseq::decode() const {
	string seq;
	seq.reserve(length());
	for(int8_t b : *this)
		seq.push_back(DNAalphabet::decode(b));
	return seq;
}

DNAseq& DNAseq::reverse() {
	std::reverse(begin(), end());
	return *this;
}

DNAseq& DNAseq::complement() {
	for(DNAseq::size_type i = 0; i < length(); ++i)
		(*this)[i] = DNAalphabet::complement((*this)[i]);
	return *this;
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
	DNAseq::reverse_iterator gap_rstart, gap_rend;
	while((gap_rstart = std::search_n(rbegin(), rend(), minNGap, DNAalphabet::N)) != rend()) { /* reverse gaps backward */
		for(gap_rend = gap_rstart; *gap_rend == DNAalphabet::N; ++gap_rend) /* find gap end on reverse order */
			continue;
		erase(gap_rend.base() + 1, gap_rstart.base());
	}
	return *this;
}

DNAseq& DNAseq::trimGaps(int mode) {
	if(mode | THREE_PRIME) { /* backward search */
		DNAseq::reverse_iterator gap_rend;
		for(gap_rend = rbegin(); DNAalphabet::isGap(*gap_rend) && gap_rend != rend(); ++gap_rend)
			continue;
		erase(gap_rend.base(), end());
	}
	if(mode | FIVE_PRIME) { /* forward search */
		DNAseq::iterator gap_end;
		for(gap_end = begin(); DNAalphabet::isGap(*gap_end) && gap_end != end(); ++gap_end)
			continue;
		erase(begin(), gap_end);
	}
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */


