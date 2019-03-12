/*
 * DNAseq.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */
#include <algorithm>
#include <cassert>
#include <htslib/sam.h>
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGseqTK {
namespace dna {

DNAseq& removeGaps(DNAseq& seq) {
	seq.erase(std::remove(seq.begin(), seq.end(), DNAalphabet::GAP_BASE), seq.end());
	return seq;
}

DNAseq encode(const string& seqStr) {
	DNAseq seq;
	seq.reserve(seqStr.length());
	for(string::value_type c : seqStr)
		seq.push_back(DNAalphabet::encode(c));
	return seq;
}

string decode(const DNAseq& seq) {
	string seqStr;
	seqStr.reserve(seq.length());
	for(DNAseq::value_type b : seq)
		seqStr.push_back(DNAalphabet::decode(b));
	return seqStr;
}

istream& read(DNAseq& seq, istream& in) {
	string seqStr;
	in >> seqStr;
	seq = encode(seqStr);
	return in;
}

DNAseq nt16Encode(const DNAseq& seq) {
	const size_t L = seq.length(); // seq length
	const size_t N = (L + 1) / 2; // compressed length
	DNAseq seqNt16(N, 0);
	for(size_t i = 0; i < L; i += 2)
		seqNt16[i / 2] = seq[i] << 4 | seq[i + 1]; // i+1 always valid since the added null at the end of *this
	return seqNt16;
}

DNAseq nt16Decode(size_t L, const DNAseq& seqNt16) {
	assert(seqNt16.length() == (L + 1) / 2);
	DNAseq seq(L, 0); // 0 init a DNAseq
	for(size_t i = 0; i < L; ++i)
		seq[i] = bam_seqi(seqNt16, i);
	return seq;
}

istream& nt16Load(DNAseq& seq, istream& in) {
	DNAseq nt16Seq;
	size_t L = 0; // original length
	in.read((char*) &L, sizeof(size_t));
	size_t N = (L + 1) / 2; // compressed length
	StringUtils::loadString(nt16Seq, in, N);
	seq = nt16Decode(L, nt16Seq);
	return in;
}

ostream& nt16Save(const DNAseq& seq, ostream& out) {
	const size_t L = seq.length(); // original length
	const size_t N = (L + 1) / 2; // compressed length
	out.write((const char*) &L, sizeof(size_t));
	return StringUtils::saveString(nt16Encode(seq), out, N);
}

BaseCount baseCount(const DNAseq& seq, size_t from, size_t len) {
	BaseCount bc;
	bc.fill(0);
	for(size_t i = from; i < from + len; ++i)
		bc[seq[i]]++;
	return bc;
}

} /* namespace dna */
} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
