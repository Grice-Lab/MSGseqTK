/*
 * RRFMIndex.cpp
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <cstddef>
#include "RRFMIndex.h"
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

RRFMIndex& RRFMIndex::build(const DNAseq& seq) throw(std::length_error) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
}

void RRFMIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::value_type b : seq)
		if(DNAalphabet::isValid(b))  /* b is in 1..5 */
			C[b]++;
}

ostream& RRFMIndex::save(ostream& out) const {
	out.write((char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bwt->save(out);
	return out;
}

istream& RRFMIndex::load(istream& in) {
	in.read((char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bwt = WaveletTreeNoptrs::load(in);
	return in;
}

saidx_t RRFMIndex::count(const DNAseq& pattern) const {
	if(pattern.empty())
		return 0; /* empty pattern matches to nothing */

    saidx_t start = 1;
    saidx_t end = bwt->getLength();
	/* search pattern left-to-right, as bwt is the reverse FM-index */
    for (DNAseq::value_type b : pattern) {
      start = C[b] + bwt->rank(b, start - 1); /* LF Mapping */
      end = C[b] + bwt->rank(b, end) - 1; /* LF Mapping */
    }
    return start <= end ? end - start + 1 : 0;
}

void RRFMIndex::buildBWT(const DNAseq& seq) {
	/* construct reverse seq */
	const size_t N =  seq.length();
	sauchar_t* rseq = new sauchar_t[N];
	std::copy(seq.rbegin(), seq.rend(), rseq);

	/* transfer rseq to bwt */
	sauchar_t* rbwt = rseq;
	divbwt(rseq, rbwt, NULL, N);

	/* construct RRR_compressed BWT */
    bwt = new WaveletTreeNoptrs((uint *) rbwt, N, sizeof(sauchar_t) * 8,
    		new BitSequenceBuilderRRR(RRR_SAMPLE_RATE), /* smart ptr */
			new MapperNone() /* smart ptr */,
			true); // free the rbwt after use
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
