/*
 * DNAseq.h
 *  a DNAseq is just a typedef with static methods to operator it
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */

#ifndef SRC_DNASEQ_H_
#define SRC_DNASEQ_H_

#include <string>
#include <iostream>
#include <algorithm>
#include <array>
#include "DNAalphabet.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::basic_string;
using std::istream;
using std::ostream;

/**
 * A DNAseq is just a convenient typedef of basic_string with supported stand-alone functions
 */
typedef std::basic_string<nt16_t> DNAseq;
typedef std::array<size_t, DNAalphabet::SIZE> BaseCount;

namespace dna {

inline
bool isValid(const DNAseq& seq, DNAseq::size_type i) {
	return DNAalphabet::isValid(seq[i]);
}

/** test whether the whole sequence is valid */
inline
bool isValid(const DNAseq& seq) {
	return std::all_of(seq.begin(), seq.end(),
			[](DNAseq::value_type b) { return DNAalphabet::isValid(b); });
}

/** test whether position i is a valid base (non-gap) */
inline
bool isBase(const DNAseq& seq, DNAseq::size_type i) {
	return DNAalphabet::isBase(seq[i]);
}

/** test whether the whole sequence is all base, no gap */
inline
bool isBase(const DNAseq& seq) {
	return std::all_of(seq.begin(), seq.end(),
			[](DNAseq::value_type b) { return DNAalphabet::isBase(b); });
}

/** test whether any base is gap */
inline
bool hasGap(const DNAseq& seq) {
	return std::any_of(seq.begin(), seq.end(),
			[](DNAseq::value_type b) { return DNAalphabet::isGap(b); });
}

/** reverse this a DNAseq */
inline
DNAseq& reverse(DNAseq& seq) {
	std::reverse(seq.begin(), seq.end());
	return seq;
}

/** Get a reverse copy of a DNAseq */
inline
DNAseq reverse(const DNAseq& seq) {
	DNAseq rSeq(seq);
	return reverse(rSeq);
}

/** Get a complement of a DNAseq */
inline
DNAseq& complement(DNAseq& seq) {
	for(DNAseq::value_type& b: seq)
		b = DNAalphabet::complement(b);
	return seq;
}

/** Get a complement copy of a DNAseq */
inline
DNAseq complement(const DNAseq& seq) {
	DNAseq cSeq(seq);
	return complement(cSeq);
}

/** reverse complement a DNAseq */
inline
DNAseq& revcom(DNAseq& seq) {
	return complement(reverse(seq));
}

/**
 * Generate a reverse complement copy of a DNAseq
 */
inline
DNAseq revcom(const DNAseq& seq) {
	DNAseq rcSeq(seq);
	return revcom(rcSeq);
}

/** transform a DNAseq to basic-only */
inline
DNAseq& toBasic(DNAseq& seq) {
	for(DNAseq::value_type& b : seq)
		b = DNAalphabet::toBasic(b);
	return seq;
}

/** get a basic-only copy of a DNAseq */
inline
DNAseq toBasic(const DNAseq& seq) {
	DNAseq bSeq(seq);
	return toBasic(bSeq);
}

DNAseq& removeGaps(DNAseq& seq);

/** encode a DNAseq from seqStr */
DNAseq encode(const string& seqStr);

/** decode a DNAseq to rawStr */
string decode(const DNAseq& seq);

/** decode a DNAseq at given position */
inline
char decode(const DNAseq& seq, size_t i) {
	return DNAalphabet::decode(seq[i]);
}

/** encode a DNAseq into an nt16 coded DNAseq */
DNAseq nt16Encode(const DNAseq& seq);

/** decode a an nt16 coded DNAseq to DNAseq */
DNAseq nt16Decode(size_t L, const DNAseq& nt16Seq);

/** load a DNAseq from binary input */
inline
istream& load(DNAseq& seq, istream& in) {
	return StringUtils::loadString(seq, in);
}

/** save a DNAseq to binary output */
inline
ostream& save(const DNAseq& seq, ostream& out) {
	return StringUtils::saveString(seq, out);
}

/** read a DNAseq from text input */
istream& read(DNAseq& seq, istream& in);

/** write a DNAseq to text output */
inline
ostream& write(const DNAseq& seq, ostream& out) {
	return out << decode(seq);
}

/** load a DNAseq from binary input with compressed nt16 encoding */
istream& nt16Load(DNAseq& seq, istream& in);

/** save a DNAseq to binary output with compressed nt16 encoding */
ostream& nt16Save(const DNAseq& seq, ostream& out);

/** get base count of a given region */
BaseCount baseCount(const DNAseq& seq, size_t start = 0, size_t len = DNAseq::npos);

/** read a DNAseq from text input */
inline
istream& operator>>(istream& in, DNAseq& seq) {
	return read(seq, in);
}

/** write a DNAseq to text output */
inline
ostream& operator<<(ostream& out, const DNAseq& seq) {
	return write(seq, out);
}

} /* namespace dna */
} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNASEQ_H_ */
