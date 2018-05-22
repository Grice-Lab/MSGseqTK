/*
 * RRFMIndex.cpp
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#include "RFMIndex.h"

#include <algorithm>
#include <cstddef>
#include "BitSequenceBuilder.h"
#include "BitSequenceBuilderRRR.h"
#include "Mapper.h"
#include "MapperNone.h"

namespace EGriceLab {
namespace MSGseqClean {

using cds_static::Mapper;
using cds_static::MapperNone;
using cds_static::BitSequenceBuilder;
using cds_static::BitSequenceBuilderRRR;
using cds_static::WaveletTreeNoptrs;
using cds_static::BitString;

RFMIndex& RFMIndex::build(const DNAseq& seq) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
}

void RFMIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::value_type b : seq)
		C[b]++;
	C['\0']++; // count null terminator
	/* construct cumulative counts (total counts smaller than this character) */
    saidx_t prev = C[0];
    saidx_t tmp;
    C[0] = 0;
    for (int i = 1; i <= DNAalphabet::SIZE; ++i) {
      tmp = C[i];
      C[i] = C[i-1] + prev;
      prev = tmp;
    }
}

ostream& RFMIndex::save(ostream& out) const {
	out.write((char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bwt->save(out);
	return out;
}

istream& RFMIndex::load(istream& in) {
	in.read((char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bwt.reset(WaveletTreeNoptrs::load(in));
	return in;
}

saidx_t RFMIndex::count(const DNAseq& pattern) const {
	if(pattern.empty())
		return 0; /* empty pattern matches to nothing */
	if(!pattern.allBase())
		return 0;

    saidx_t start = 1; /* 1-based start */
    saidx_t end = length() - 1; /* 1-based end */
	/* search pattern left-to-right, as bwt is the reverse FM-index */
    for(DNAseq::const_reverse_iterator b = pattern.rbegin(); b != pattern.rend() && start <= end; ++b) {
    	start = LF(*b, start - 1); /* LF Mapping */
    	end = LF(*b, end) - 1; /* LF Mapping */
    }
    return start <= end ? end - start + 1 : 0;
}

RFMIndex& RFMIndex::operator+=(const RFMIndex& other) {
	if(bwt == other.bwt) /* prevent self addition */
		return *this;
	if(!isInitiated())
		return *this = other;

	/* build interleaving bitvector between *this and other */
	const saidx_t N1 = length();
	const saidx_t N2 = other.length();
	const saidx_t N = N1 + N2;

	/* build merged C[] */
	saidx_t CMerged[UINT8_MAX + 1] = { 0 };
	for(uint8_t i = 0; i < UINT8_MAX; ++i)
		CMerged[i] = C[i] + other.C[i];

	/* build RA */
	saidx_t* RA = new saidx_t[N1];

	for(saidx_t i = 0, j = 0; (j = LF(i) - 1) != 0; i = j) {
		saidx_t val = other.LF(bwt->access(i), val) - 1;
		RA[j] = val;
		i = j;
	}

	for(saidx_t i = 0; i < N1; ++i)
		cerr << "RA[" << i << "]: " << RA[i] << endl;

	/* build interleaving bitvector */
	BitString B(N);
	for(saidx_t i = 0; i < N1; ++i)
		B.setBit(i + 1 + RA[i]);

	cerr << "B: ";
	for(saidx_t i = 0; i < N; ++i)
		cerr << B.getBit(i) << " ";
	cerr << endl;

	/* build merbed BWT */
	sauchar_t* bwtNew = new sauchar_t[N];
	for(saidx_t i = 0, j = 0, k = 0; k < N; ++k) {
		bwtNew[k] = B.getBit(k) /* use this or other */ ? bwt->access(i++) : other.bwt->access(j++);
	}
    BWTRRR_ptr bwtMerged = std::make_shared<BWTRRR>(reinterpret_cast<uint*> (bwtNew), N, sizeof(sauchar_t) * 8,
    		new BitSequenceBuilderRRR(RRR_SAMPLE_RATE), /* smart ptr */
			new MapperNone() /* smart ptr */,
			false); // do not free the bwtNew after use
    delete[] bwtNew;

    /* swap data */
    std::swap(C, CMerged);
    std::swap(bwt, bwtMerged);

	return *this;
}

string RFMIndex::getBWT() const {
	string bwtStr;
	bwtStr.reserve(length());
	for(saidx_t i = 0; i < length(); ++i) {
		uint b = bwt->access(i);
		bwtStr.push_back(b != '\0' ? DNAalphabet::decode(b) : TERMINAL_SYMBOL);
	}
	return bwtStr;
}

DNAseq RFMIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	DNAseq seq;
	seq.reserve(length() - 1);
	for(saidx_t i = 0; seq.length() < length() - 1; i = LF(i) - 1)
		seq.push_back(bwt->access(i));
	std::reverse(seq.begin(), seq.end());
	return seq;
}

void RFMIndex::buildBWT(const DNAseq& seq) {
	const size_t N =  seq.length() + 1;
	/* construct SA */
    saidx_t* SA = new saidx_t[N];
    saidx_t errn = divsufsort((const sauchar_t*) seq.c_str(), SA, N);
	if(errn != 0)
		throw runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build bwt */
	sauchar_t* bwtSeq = new sauchar_t[N];
    for(saidx_t i = 0; i < N; ++i)
        if(SA[i] == 0) // matches to the null
            bwtSeq[i] = 0; // null terminal
        else bwtSeq[i] = seq[SA[i] - 1];
//	divbwt((const sauchar_t*) seq.c_str(), bwtSeq, SA, N);

	/* construct RRR_compressed BWT */
    bwt = std::make_shared<BWTRRR>(reinterpret_cast<uint*> (bwtSeq), N, sizeof(sauchar_t) * 8,
    		new BitSequenceBuilderRRR(RRR_SAMPLE_RATE), /* smart ptr */
			new MapperNone() /* smart ptr */,
			false); // do not free the rbwt after use

    delete[] bwtSeq;
    delete[] SA;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
