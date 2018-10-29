/*
 * RRFMIndex.cpp
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <utility>
#include <cstddef>
#include <cassert>
#include <stdexcept>
#include "FMIndex.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::vector;
using EGriceLab::libSDS::BitStr32;

const saidx_t FMIndex::totalBases() const {
	saidx_t N = 0;
	for(int8_t i = DNAalphabet::A; i < DNAalphabet::SIZE; i++)
		N += B[i];
	return N;
}

FMIndex& FMIndex::build(const DNAseq& seq) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
}

void FMIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::value_type b : seq)
		B[b]++;
	B[0]++; // count terminal null

	/* calculate cumulative counts (total counts smaller than this character) */
    saidx_t sum = 0;
    for(int i = 0; i <= DNAalphabet::SIZE; ++i) {
    	C[i] = sum; /* use sum before adding current B[i] */
    	sum += B[i];
    }
}

ostream& FMIndex::save(ostream& out) const {
	out.write((const char*) B, (UINT8_MAX + 1) * sizeof(saidx_t));
	out.write((const char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bwt.save(out);
	size_t nSAsampled = SAsampled.size();
	out.write((const char*) &nSAsampled, sizeof(size_t));
	if(nSAsampled) {
		assert(nSAsampled == length() / SA_SAMPLE_RATE + C[0 + 1] + 1);
		out.write((const char*) SAsampled.data(), sizeof(saidx_t) * nSAsampled);
		SAbit.save(out);
	}
	return out;
}

istream& FMIndex::load(istream& in) {
	in.read((char*) B, (UINT8_MAX + 1) * sizeof(saidx_t));
	in.read((char*) C, (UINT8_MAX + 1) * sizeof(saidx_t));
	bwt.load(in);
	size_t nSAsampled = 0;
	in.read((char*) &nSAsampled, sizeof(size_t));
	if(nSAsampled) {
		assert(nSAsampled == length() / SA_SAMPLE_RATE + C[0 + 1] + 1);
		SAsampled.resize(nSAsampled);
		in.read((char*) SAsampled.data(), sizeof(saidx_t) * nSAsampled);
		SAbit.load(in);
	}
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
		*this = other;
		return *this;
	}

	/* build interleaving bitvector between other and *this */
	const saidx_t N1 = other.length();
	const saidx_t N2 = length();
	const saidx_t N = N1 + N2;

	/* build RA and interleaving bitvector */
	BitStr32 bstr(N);
	for(saidx_t i = 0, j = 0, shift = 0, RA = C[0 + 1]; j < N1; ++j) {
		bstr.set(i + RA);
		/* LF mapping */
		sauchar_t b = other.bwt.access(i);
		if(b == 0) {
			RA = C[0 + 1];
			i = ++shift;
		}
		else {
			RA = LF(b, RA - 1);
			i = other.LF(i) - 1;
		}
//		cerr << "i: " << i << " c: " << DNAalphabet::decode(b) << " RA: " << RA << endl;
	}
	std::cerr << "bstr constructed" << std::endl;

	/* build merbed BWT */
	basic_string<sauchar_t> bwtM;
	bwtM.reserve(N);
	for(saidx_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM.push_back(bstr.test(k) ? other.bwt.access(i++) : bwt.access(j++));
	std::cerr << "bwtM constructed" << std::endl;
	/* update bwtRRR */
    bwt = WaveletTreeRRR(bwtM, DNAalphabet::N, DNAalphabet::T);
	std::cerr << "bwt constructed" << std::endl;
    bwtM.clear();

	/* merging B[] and C[] */
	for(saidx_t i = 0; i <= DNAalphabet::SIZE; ++i) {
		B[i] += other.B[i];
		C[i] += other.C[i];
	}
	std::cerr << "B and C merged" << std::endl;

    /* reset SA */
	SAsampled.clear();
	SAbit.reset();

    if(keepSA)
    	buildSA();

	return *this;
}

DNAseq FMIndex::getBWT() const {
	DNAseq bwtSeq;
	bwtSeq.reserve(length());
	for(saidx_t i = 0; i < length(); ++i)
		bwtSeq.push_back(bwt.access(i));
	return bwtSeq;
}

DNAseq FMIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	DNAseq seq;
	seq.reserve(length() - 1);
	for(saidx_t i = 0, shift = 0; seq.length() < length() - 1;) {
		sauchar_t b = bwt.access(i);
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
		throw std::runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build bwt */
	basic_string<sauchar_t> bwtSeq;
	bwtSeq.reserve(N);

    for(saidx_t i = 0; i < N; ++i)
    	bwtSeq.push_back(SA[i] == 0 ? 0 : seq[SA[i] - 1]);
    delete[] SA;

	/* construct BWTRRR */
    bwt = WaveletTreeRRR(bwtSeq, DNAalphabet::N, DNAalphabet::T);
    bwtSeq.clear();

	if(keepSA) /* if intermediate SA need to be kept */
		buildSA();
}

void FMIndex::buildSA() {
	assert(isInitiated());
	assert(SAsampled.empty());
	const saidx_t N = length();
	/* build a BitStr in the 1st pass */
	BitStr32 bstr(N);
	for(saidx_t i = 0; i < N; ++i) {
		if(bwt.access(i) == 0 || i % SA_SAMPLE_RATE == 0)
			bstr.set(i);
	}
	std::cerr << "bstr built" << std::endl;
	SAbit = BitSeqGGMN(bstr); /* reset the SAbit */
	std::cerr << "SAbit built" << std::endl;

	/* build SAsampled in the 2nd pass */
	std::cerr << "SAsampled.size(): " << SAsampled.size() << std::endl;
	std::cerr << "SAsampled resizing to " << (N / SA_SAMPLE_RATE + C[0 + 1] + 1) << std::endl;
	SAsampled.reserve(N / SA_SAMPLE_RATE + C[0 + 1] + 1);
	SAsampled.resize(N / SA_SAMPLE_RATE + C[0 + 1] + 1); // plus all Ns */
	std::cerr << "SAsampled resized to " << SAsampled.size() << std::endl;
	for(saidx_t i = 0, j = N, shift = 0; j > 0; --j) { /* i: 0-based on BWT; j: 1-based on seq */
		sauchar_t b = bwt.access(i);
		if(b == 0 || i % SA_SAMPLE_RATE == 0)
			SAsampled[SAbit.rank1(i) - 1] = j - 1;
		/* LF-mapping */
		i = b != 0 ? LF(i) - 1 : ++shift;
	}
	std::cerr << "SAsampled built" << std::endl;
}

saidx_t FMIndex::accessSA(saidx_t i) const {
	saidx_t dist = 0;
	while(!SAbit.access(i)) {
		i =  LF(i) - 1; // backward LF-mapping
		dist++;
	}
	return SAsampled[SAbit.rank1(i) - 1] + dist;
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


} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
