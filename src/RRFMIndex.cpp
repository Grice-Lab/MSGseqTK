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

RRFMIndex& RRFMIndex::build(const DNAseq& seq) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
}

void RRFMIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::value_type b : seq)
		C[b]++;
	C[DNAalphabet::N]++; /* add a terminal gap */

	/* construct cumulative counts, excluding the gap char */
    saidx_t prev = C[0];
    saidx_t tmp;
    C[0] = 0;
    for (int i = 1; i < DNAalphabet::SIZE; ++i) {
      tmp = C[i];
      C[i] = C[i-1] + prev;
      prev = tmp;
    }
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
	if(!pattern.allBase())
		return 0;

    saidx_t start = 1;
    saidx_t end = bwt->getLength();
	/* search pattern left-to-right, as bwt is the reverse FM-index */
    for (DNAseq::const_iterator b = pattern.begin(); b != pattern.end() && start <= end; ++b) {
    	start = C[*b] + bwt->rank(*b, start - 1); /* LF Mapping */
    	end = C[*b] + bwt->rank(*b, end) - 1; /* LF Mapping */
    }
    return start <= end ? end - start + 1 : 0;
}

void RRFMIndex::buildBWT(const DNAseq& seq) {
	/* construct reverse seq */
	const size_t N =  seq.length() + 1; /* null/N terminated */
	sauchar_t* rseq = new sauchar_t[N]();
	std::copy(seq.rbegin(), seq.rend(), rseq);
	rseq[N - 1] = 0; /* null terminator */

	/* construct SA */
    saidx_t* SA = new saidx_t[N];
    saidx_t errn = divsufsort(rseq, SA, N);
	if(errn != 0)
		throw runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build bwt */
	sauchar_t* rbwt = rseq; /* transfer in-place */
	divbwt(rseq, rbwt, SA, N);

	/* construct RRR_compressed BWT */
    bwt = new WaveletTreeNoptrs(reinterpret_cast<uint*> (rbwt), N, sizeof(sauchar_t) * 8,
    		new BitSequenceBuilderRRR(RRR_SAMPLE_RATE), /* smart ptr */
			new MapperNone() /* smart ptr */,
			false); // do not free the rbwt after use

    delete[] rseq;
    delete[] SA;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
