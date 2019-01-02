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
const DNAseq DNAseq::DNAgap(1, DNAalphabet::GAP_BASE); /* single base gap DNAseq */

string DNAseq::decode() const {
	string seq;
	seq.reserve(length());
	for(DNAseq::value_type b : *this)
		seq.push_back(DNAalphabet::decode(b));
	return seq;
}

DNAseq& DNAseq::reverse() {
	std::reverse(begin(), end());
	return *this;
}

DNAseq& DNAseq::complement() {
	for(DNAseq::value_type& b: *this)
		b = DNAalphabet::complement(b);
	return *this;
}

DNAseq& DNAseq::assign(const string& str) {
	resize(str.length());
	std::transform(str.begin(), str.end(), begin(), DNAalphabet::encode);
	return *this;
}

DNAseq& DNAseq::append(const string& str) {
	reserve(length() + str.length());
	for(char s : str)
		push_back(DNAalphabet::encode(s));
	return *this;
}

istream& DNAseq::read(istream& in) {
	string str;
	in >> str;
	assign(str);
	return in;
}

DNAseq& DNAseq::removeInvalid() {
	erase(std::remove_if(begin(), end(),
			[](DNAseq::value_type b) { return !DNAalphabet::isValid(b); }),
			end());
	return *this;
}

DNAseq& DNAseq::removeGaps() {
	erase(std::remove(begin(), end(), DNAalphabet::GAP_BASE), end());
	return *this;
}

DNAseq DNAseq::substr(size_t pos, size_t len) const {
	DNAseq seg;
	if(pos + len >= length())
		len = length() - pos;
	seg.resize(len);
	std::copy_n(begin() + pos, len, seg.begin());
	return seg;
}

BAM::seq_str DNAseq::nt16Encode() const {
	const size_t L = length(); // raw length
	const size_t L16 = (L + 1) / 2; // ceil(L / 2)
	BAM::seq_str seqNt16(L16, 0);
	for(size_t i = 0; i < L; i += 2)
		seqNt16[i / 2] = (*this)[i] << 4 | (*this)[i + 1]; // i+1 always valid since the added null at the end of *this
	return seqNt16;
}

DNAseq DNAseq::nt16Decode(size_t L, const BAM::seq_str& seqNt16) {
	assert(seqNt16.length() == (L + 1) / 2);
	DNAseq seq(L, 0); // 0 init a DNAseq
	for(size_t i = 0; i < L; ++i)
		seq[i] = BAM::getSeqBase(seqNt16, i);
	return seq;
}

istream& DNAseq::nt16Load(istream& in) {
	size_t L = 0; // uncompressed length
	in.read((char*) &L, sizeof(size_t));
	resize(L);
	for(size_t i = 0; i < L; i += 2) {
		value_type b = in.get();
		(*this)[i] = (b & DNAalphabet::NT16_UPPER_MASK) >> 4;
		(*this)[i + 1] = b & DNAalphabet::NT16_LOWER_MASK; // i+1 always valid
	}
	return in;
}

ostream& DNAseq::nt16Save(ostream& out) const {
	const size_t L = length(); // raw length
	out.write((const char*) &L, sizeof(size_t));
	for(size_t i = 0; i < L; i += 2) { // (*this)[i+1] always valid
		out.put((*this)[i] << 4 | (*this)[i + 1]); // value_type can be fit into char
	}
	return out;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
