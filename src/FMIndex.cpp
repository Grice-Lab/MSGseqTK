/*
 * RRFMIndex.cpp
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#include "FMIndex.h"

#include <algorithm>
#include <cstddef>
#include <cassert>
#include "BitSequenceBuilder.h"
#include "BitSequenceBuilderRRR.h"
#include "Mapper.h"
#include "MapperNone.h"

namespace EGriceLab {
namespace MSGseqClean {
using std::vector;
using cds_static::Mapper;
using cds_static::MapperNone;
using cds_static::BitSequenceRRR;
using cds_static::BitSequenceBuilder;
using cds_static::BitSequenceBuilderRRR;
using cds_static::WaveletTreeNoptrs;
using cds_static::BitString;

FMIndex& FMIndex::build(const DNAseq& seq) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
}

void FMIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::const_iterator b = seq.begin(); b < seq.end(); ++b)
		C[*b]++;
	C[0]++; // always terminated with null
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

ostream& FMIndex::save(ostream& out) const {
	out.write((const char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bool SAflag = !SAsampled.empty() && SAbit != nullptr;
	out.write((const char*) &SAflag, sizeof(bool));
	if(SAflag) {
		assert(SAsampled.size() == length() / SA_SAMPLE_RATE + C[0 + 1] + 1);
		StringUtils::saveString(SAsampled, out);
		SAbit->save(out);
	}
	bwt->save(out);
	return out;
}

istream& FMIndex::load(istream& in) {
	in.read((char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bool SAflag;
	in.read((char*) &SAflag, sizeof(bool));
	if(SAflag) {
		StringUtils::loadString(SAsampled, in);
		SAbit.reset(BitSequenceRRR::load(in));
	}

	bwt.reset(WaveletTreeNoptrs::load(in));
	return in;
}

saidx_t FMIndex::count(const DNAseq& pattern) const {
	if(pattern.empty())
		return 0; /* empty pattern matches to nothing */
	if(!pattern.allBase())
		return 0;

    saidx_t start = 0;
    saidx_t end = length() - 1;
	/* search pattern left-to-right, as bwt is the reverse FM-index */
    for(DNAseq::const_reverse_iterator b = pattern.rbegin(); b != pattern.rend() && start <= end; ++b) {
    	if(start == 0) {
    		start = C[*b];
    		end = C[*b + 1] - 1;
    	}
    	else {
    		start = LF(*b, start - 1); /* LF Mapping */
    		end = LF(*b, end) - 1; /* LF Mapping */
    	}
    }
    return start <= end ? end - start + 1 : 0;
}

FMIndex& FMIndex::operator+=(const FMIndex& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other; /* shallow copy the object */
		return *this;
	}

	/* build interleaving bitvector between *this and other */
	const saidx_t N1 = length();
	const saidx_t N2 = other.length();
	const saidx_t N = N1 + N2;

	/* build merged C[] */
	saidx_t CMerged[UINT8_MAX + 1] = { 0 };
	for(saidx_t i = 0; i <= DNAalphabet::SIZE; ++i)
		CMerged[i] = C[i] + other.C[i];

	/* build RA and interleaving bitvector */
	BitString B(N);
	for(saidx_t i = 0, j = 0, shift = 0, RA = other.C[0 + 1]; j < N1; ++j) {
		B.setBit(i + RA);
		/* LF mapping */
		sauchar_t b = bwt->access(i);
		if(b == 0) {
			RA = other.C[0 + 1];
			i = ++shift;
		}
		else {
			RA = other.LF(b, RA - 1);
			i = LF(i) - 1;
		}
//		cerr << "i: " << i << " c: " << DNAalphabet::decode(b) << " RA: " << RA << endl;
	}

	/* build merbed BWT */
	sauchar_t* bwtM= new sauchar_t[N];
	for(saidx_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM[k] = B.getBit(k) ? bwt->access(i++) : other.bwt->access(j++);

    BWTRRR_ptr bwtMerged = std::make_shared<BWTRRR>(reinterpret_cast<uint*> (bwtM), N, sizeof(sauchar_t) * 8,
    		new BitSequenceBuilderRRR(RRR_SAMPLE_RATE), /* smart ptr */
			new MapperNone() /* smart ptr */,
			false); // do not free the bwtNew after use
    delete[] bwtM;

    /* swap data */
    std::swap(C, CMerged);
    std::swap(bwt, bwtMerged);

    if(keepSA)
    	buildSA();

	return *this;
}

DNAseq FMIndex::getBWT() const {
	DNAseq bwt;
	bwt.reserve(length());
	for(saidx_t i = 0; i < length(); ++i)
		bwt.push_back(this->bwt->access(i));
	return bwt;
}

DNAseq FMIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	DNAseq seq;
	seq.reserve(length() - 1);
	for(saidx_t i = 0, shift = 0; seq.length() < length() - 1;) {
		sauchar_t b = bwt->access(i);
		seq.push_back(b);
		i = b != 0 ? LF(i) - 1 : ++shift;
	}
	std::reverse(seq.begin(), seq.end());
	return seq;
}

void FMIndex::buildBWT(const DNAseq& seq) {
	const size_t N =  seq.length() + 1;
	/* construct SA */
    saidx_t* SA = new saidx_t[N];
    saidx_t errn = divsufsort((const sauchar_t*) (seq.c_str()), SA, N);
	if(errn != 0)
		throw runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build bwt */
	sauchar_t* bwtSeq = new sauchar_t[N];
#pragma omp parallel for
    for(saidx_t i = 0; i < N; ++i) {
        if(SA[i] == 0) // matches to the null
            bwtSeq[i] = 0; // null terminal
        else bwtSeq[i] = seq[SA[i] - 1];
    }

	/* construct RRR_compressed BWT */
    bwt = std::make_shared<BWTRRR>(reinterpret_cast<uint*> (bwtSeq), N, sizeof(sauchar_t) * 8,
    		new BitSequenceBuilderRRR(RRR_SAMPLE_RATE), /* smart ptr */
			new MapperNone() /* smart ptr */,
			false); // do not free the rbwt after use

    delete[] bwtSeq;
    delete[] SA;

	if(keepSA) /* if intermediate SA need to be kept */
		buildSA();
}

void FMIndex::buildSA() {
	assert(isInitiated());
	const saidx_t N = length();
	/* build BitVector in the 1st pass */
	{
		BitString B(N);
		for(saidx_t i = 0; i < N; ++i) {
			if(bwt->access(i) == 0 || i % SA_SAMPLE_RATE == 0)
				B.setBit(i);
		}
		SAbit.reset(new BitSequenceRRR(B, RRR_SAMPLE_RATE)); /* use RRR implementation */
	}

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(N / SA_SAMPLE_RATE + C[0 + 1] + 1); /* need additional storage for all Ns */
	for(saidx_t i = 0, j = N, shift = 0; j > 0; --j) { /* i: 0-based on BWT; j: 1-based on seq */
		sauchar_t b = bwt->access(i);
		if(b == 0 || i % SA_SAMPLE_RATE == 0)
			SAsampled[SAbit->rank1(i) - 1] = j - 1;
		/* LF-mapping */
		i = b != 0 ? LF(i) - 1 : ++shift;
	}
}

saidx_t FMIndex::accessSA(saidx_t i) const {
	saidx_t dist = 0;
	while(!SAbit->access(i)) {
		i =  LF(i) - 1; // backward LF-mapping
		dist++;
	}
	return SAsampled[SAbit->rank1(i) - 1] + dist;
}

vector<Loc> FMIndex::locateAll(const DNAseq& pattern) const {
	vector<Loc> locs;
	if(pattern.empty())
		return locs;
	if(!pattern.allBase())
		return locs;

    saidx_t start = 0; /* 0-based */
    saidx_t end = length() - 1; /* 1-based */
	/* backward pattern search */
    for(DNAseq::const_reverse_iterator b = pattern.rbegin(); b != pattern.rend() && start <= end; ++b) {
    	if(start == 0) {
    		start = C[*b];
    		end = C[*b + 1] - 1;
    	}
    	else {
    		start = LF(*b, start - 1); /* LF Mapping */
    		end = LF(*b, end) - 1; /* LF Mapping */
    	}
    }

    for(saidx_t i = start; i <= end; ++i) {
    	saidx_t SAstart = accessSA(i);
    	locs.push_back(Loc(SAstart, SAstart + pattern.length()));
    }
    return locs;
}

Loc FMIndex::reverseLoc(const Loc& loc) const {
	return Loc(reverseLoc(loc.end - 1) - 1, reverseLoc(loc.start));
}

MEM FMIndex::getMEM(const DNAseq&read, saidx_t from) const {
	uint64_t start = 0;
	uint64_t end = 0;
	uint64_t nextStart = start;
	uint64_t nextEnd = end;
	uint64_t to;
	/* search read left-to-right */
	for(to = from; to < read.length(); ++to, start = nextStart, end = nextEnd) {
		sauchar_t b = read[to];
		if(b == 0) /* null gap */
			break;
		if(start == 0) {
			nextStart = C[b];
			nextEnd = C[b + 1] - 1;
		}
		else {
			nextStart = LF(b, start - 1);
			nextEnd = LF(b, end) - 1;
		}

		if(nextStart > nextEnd)
			break;
	}

	/* construct MEM with basic matching info */
	MEM mem(from, to, &read, nullptr);
	if(start == 0 && end == 0)
		return mem;
	/* adding all matching locs */
	for(saidx_t i = start; i <= end; ++i) {
		saidx_t SAstart = accessSA(i);
		mem.locs.push_back(Loc(SAstart, SAstart + to - from));
	}

	return mem;
}

MEM FMIndex::getMEM(const DNAseq&read, const QualStr& qual, saidx_t from) const {
	MEM mem = getMEM(read, from);
	mem.qual = &qual;
	return mem;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
