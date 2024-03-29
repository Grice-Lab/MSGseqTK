/*
 * FMDIndex.cpp
 *
 *  Created on: Apr 26, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <utility>
#include <cstddef>
#include <cassert>
#include <stdexcept>
#include "StringUtils.h"
#include "FMDIndex.h"
#include "BitSeqGGMN.h"
#include "MSGseqTK_main.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using EGriceLab::libSDS::BitStr32;
using EGriceLab::libSDS::BitSeqGGMN; // useful for temporary and fast BitSeq solution

FMDIndex::FMDIndex(const DNAseq& seq, bool buildSAsampled, int saSampleRate) {
	const int64_t N = seq.length();
	if(N > INT32_MAX) { /* use 64bit build */
		/* build SA */
		int64_t* SA = new int64_t[N];
		saint_t errn = divsufsort64((const sauchar_t*) seq.c_str(), (saidx64_t*) SA, (saidx64_t) N); // string.c_str guarantee a null terminal
		if(errn != 0)
			throw std::runtime_error("Error: Cannot build 64 bit suffix-array on input DNAseq");
		/* build counts */
		buildCounts(seq);
		/* build BWT */
		DNAseq bwt = buildBWT(seq, SA);
		/* build Gap */
		buildGap(SA);
		if(buildSAsampled)
			sampleSA(SA, bwt, saSampleRate);
		delete[] SA;
		/* build BWTRRR */
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);
	}
	else { /* use 32bit build */
		/* build SA */
		int32_t* SA = new int32_t[N];
		saint_t errn = divsufsort((const sauchar_t*) seq.c_str(), (saidx_t*) SA, (saidx_t) N); // string.c_str guarantee a null terminal
		if(errn != 0)
			throw std::runtime_error("Error: Cannot build 32 bit suffix-array on input DNAseq");
		/* build counts */
		buildCounts(seq);
		/* build BWT */
		DNAseq bwt = buildBWT(seq, SA);
		/* build Gap */
		buildGap(SA);
		/* sample SA */
		if(buildSAsampled)
			sampleSA(SA, bwt, saSampleRate);
		delete[] SA;
		/* build BWTRRR */
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);
	}
}

FMDIndex::FMDIndex(DNAseq&& seq, bool buildSAsampled, int saSampleRate) {
	const int64_t N = seq.length();
	if(N > INT32_MAX) { /* use 64bit build */
		/* build SA */
		int64_t* SA = new int64_t[N];
		saint_t errn = divsufsort64((const sauchar_t*) seq.c_str(), (saidx64_t*) SA, N); // string.c_str guarantee a null terminal
		if(errn != 0)
			throw std::runtime_error("Error: Cannot build 64 bit suffix-array on input DNAseq");
		/* build counts */
		buildCounts(seq);
		/* build BWT */
		DNAseq bwt = buildBWT(seq, SA);
		/* clear seq */
		seq.clear();
		seq.shrink_to_fit();
		/* build Gap */
		buildGap(SA);
		/* sample SA */
		if(buildSAsampled)
			sampleSA(SA, bwt, saSampleRate);
		delete[] SA;
		/* build BWTRRR */
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);
	}
	else { /* use 32bit build */
		/* build SA */
		int32_t* SA = new int32_t[N];
		saint_t errn = divsufsort((const sauchar_t*) seq.c_str(), (saidx_t*) SA, N); // string.c_str guarantee a null terminal
		if(errn != 0)
			throw std::runtime_error("Error: Cannot build 32 bit suffix-array on input DNAseq");
		/* build counts */
		buildCounts(seq);
		/* build BWT */
		DNAseq bwt = buildBWT(seq, SA);
		/* clear seq */
		seq.clear();
		seq.shrink_to_fit();
		/* build Gap */
		buildGap(SA);
		/* sample SA */
		if(buildSAsampled)
			sampleSA(SA, bwt, saSampleRate);
		delete[] SA;
		/* build BWTRRR */
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);
	}
}

void FMDIndex::buildCounts(const DNAseq& seq) {
	assert(seq.back() == 0);
	const int64_t N = seq.length();
	for(nt16_t b : seq) // one pass seq always null
		B[b]++;

	/* calculate cumulative counts */
    int64_t S = 0;
    for(nt16_t i = 0; i <= DNAalphabet::SIZE; ++i) {
    	C[i] = S;
    	S += B[i];
    }
    if(!isBiDirectional())
    	throw std::invalid_argument("input seq is not bi-directional");
    if(getExtBaseCount() > 0)
    	throw std::invalid_argument("input seq does not allow IUPAC-extended bases");
}

bool FMDIndex::isBiDirectional() const {
	for(nt16_t b = 0; b < DNAalphabet::SIZE; ++b)
		if(getBaseCount(b) != getBaseCount(DNAalphabet::complement(b)))
				return false;
	return true;
}

ostream& FMDIndex::save(ostream& out) const {
	out.write((const char*) B.data(), B.size() * sizeof(BCarray_t::value_type));
	out.write((const char*) C.data(), C.size() * sizeof(BCarray_t::value_type));
	bwtRRR.save(out);
	StringUtils::saveString(gapSA, out);
	SAidx.save(out);
	StringUtils::saveString(SAsampled, out);
	return out;
}

istream& FMDIndex::load(istream& in) {
	in.read((char*) B.data(), B.size() * sizeof(BCarray_t::value_type));
	in.read((char*) C.data(), C.size() * sizeof(BCarray_t::value_type));
	bwtRRR.load(in);
	StringUtils::loadString(gapSA, in);
	SAidx.load(in);
	StringUtils::loadString(SAsampled, in);
	return in;
}

int64_t FMDIndex::count(const DNAseq& pattern) const {
	const size_t L = pattern.length();
	size_t i = L;
	nt16_t b = pattern[i - 1];
	int64_t p = C[b];
	int64_t q = C[DNAalphabet::complement(b)];
	int64_t s = C[b + 1] - C[b];

	/* backward search */
    for(i = L - 1; i > 0 && backExt(p, q, s, pattern[i - 1]) > 0; --i)
    	continue;

    return std::max<int64_t>(s, 0);
}

FMDIndex& FMDIndex::append(const FMDIndex& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

    /* merge gap info */
    gapSA = mergeGap(*this, other);

	/* merge BWTs */
#ifndef _OPENMP
	{
		DNAseq bwt = mergeBWT(*this, other, buildInterleavingBS(*this, other));
		bwtRRR.reset();
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#else
	int nThreads = 1;
#pragma omp parallel
	nThreads = omp_get_num_threads();
	if(1 == nThreads) { // no parallelzation needed
		DNAseq bwt = mergeBWT(*this, other, buildInterleavingBS(*this, other));
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
	else {
		DNAseq bwt = mergeBWT(*this, other, BitSeqGGMN(buildInterleavingBS(*this, other)));
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#endif

	/* merge counts */
    mergeCount(other);

    /* clear SA */
    clearSA();
	return *this;
}

FMDIndex& FMDIndex::append(FMDIndex&& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

    /* merge gap info */
    gapSA = mergeGap(*this, other);

	/* merge BWTs */
#ifndef _OPENMP
	{
		DNAseq bwt = mergeBWT(*this, other, buildInterleavingBS(*this, other));
		/* reset both this and othre bwtRRR */
		other.bwtRRR.reset();
		bwtRRR.reset();
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#else
	int nThreads = 1;
#pragma omp parallel
	nThreads = omp_get_num_threads();
	if(1 == nThreads) { // no parallelzation needed
		DNAseq bwt = mergeBWT(*this, other, buildInterleavingBS(*this, other));
		other.bwtRRR.reset();
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
	else {
		DNAseq bwt = mergeBWT(*this, other, BitSeqGGMN(buildInterleavingBS(*this, other)));
		other.bwtRRR.reset();
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#endif

	/* merge counts */
    mergeCount(other);

    /* clear SA */
    clearSA();

	return *this;
}

FMDIndex& FMDIndex::prepend(const FMDIndex& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

    /* merge gap info */
    gapSA = mergeGap(other, *this);

	/* merge BWTs */
#ifndef _OPENMP
	{
		DNAseq bwt = mergeBWT(other, *this, buildInterleavingBS(other, *this));
		bwtRRR.reset();
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#else
	int nThreads = 1;
#pragma omp parallel
	nThreads = omp_get_num_threads();
	if(1 == nThreads) { // no parallelzation needed
		DNAseq bwt = mergeBWT(other, *this, buildInterleavingBS(other, *this));
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
	else {
		DNAseq bwt = mergeBWT(other, *this, BitSeqGGMN(buildInterleavingBS(other, *this)));
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#endif

	/* merge counts */
    mergeCount(other);

    /* clear SA */
    clearSA();

	return *this;
}

FMDIndex& FMDIndex::prepend(FMDIndex&& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

    /* merge gap info */
    gapSA = mergeGap(other, *this);

	/* merge BWTs */
#ifndef _OPENMP
	{
		DNAseq bwt = mergeBWT(other, *this, buildInterleavingBS(other, *this));
		other.bwtRRR.reset();
		bwtRRR.reset();
		bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#else
	int nThreads = 1;
#pragma omp parallel
	nThreads = omp_get_num_threads();
	if(1 == nThreads) { // no parallelzation needed
		DNAseq bwt = mergeBWT(other, *this, buildInterleavingBS(other, *this));
		other.bwtRRR.reset();
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
	else {
		DNAseq bwt = mergeBWT(other, *this, BitSeqGGMN(buildInterleavingBS(other, *this)));
		other.bwtRRR.reset();
		bwtRRR.reset();
	    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR
	}
#endif

	/* merge counts */
    mergeCount(other);

    /* clear SA */
    clearSA();

	return *this;
}

DNAseq FMDIndex::getBWT() const {
	const int64_t N = length();
	DNAseq bwtSeq(N, 0); // 0-init
#pragma omp parallel for
	for(int64_t i = 0; i < N; ++i)
		bwtSeq[i] = accessBWT(i);
	return bwtSeq;
}

DNAseq FMDIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	const int64_t N = length();
	DNAseq seq(N, 0);
#pragma omp parallel for schedule(dynamic)
	for(int64_t i = 0; i < numGaps(); ++i) { // i-th BWT segment
		int64_t sa = getGapSA(i); // end SA value
		seq[sa] = 0;
		nt16_t b = accessBWT(i);
		sa--;
		for(int64_t j = i; b != 0; --sa) {
			seq[sa] = b;
			j = LF(b, j) - 1; // LF-mapping
			b = accessBWT(j); // update b
		}
	}
//	assert(seq.length() == N);
	return seq;
}

DNAseq FMDIndex::buildBWT(const DNAseq& seq, const int64_t* SA) {
	assert(seq.back() == 0); // must be null-terminated
	const int64_t N = seq.length();
	/* build uncompressed bwt */
	DNAseq bwt(N, 0);
	for(size_t i = 0; i < N; ++i)
		bwt[i] = SA[i] == 0 ? 0 : seq[SA[i] - 1];
	return bwt;
}

DNAseq FMDIndex::buildBWT(const DNAseq& seq, const int32_t* SA) {
	assert(seq.back() == 0); // must be null-terminated
	const int64_t N = seq.length();
	assert(N <= INT32_MAX);
	/* build uncompressed bwt */
	DNAseq bwt(N, 0);
	for(size_t i = 0; i < N; ++i)
		bwt[i] = SA[i] == 0 ? 0 : seq[SA[i] - 1];
	return bwt;
}

void FMDIndex::buildGap(const int64_t* SA) {
	gapSA.resize(numGaps());
	for(int64_t i = 0; i < numGaps(); ++i)
		gapSA[i] = SA[i];
}

void FMDIndex::buildGap(const int32_t* SA) {
	gapSA.resize(numGaps());
	for(int64_t i = 0; i < numGaps(); ++i)
		gapSA[i] = SA[i];
}

FMDIndex& FMDIndex::sampleSA(const int64_t* SA, int64_t N, int saSampleRate) {
	BitStr32 bstr(N);
	SAsampled.clear();
	SAsampled.reserve((N + 1) / saSampleRate /* ceil(N / saSampleRate) */ + numGaps()); // enough to hold both peroid sampling and gaps
	for(size_t i = 0, j = 0; i < N; ++i) {
		if(i % saSampleRate == 0 || accessBWT(i) == 0) { /* sample at all null characters */
			bstr.set(i);
			SAsampled[j++] = SA[i];
		}
	}
	SAsampled.shrink_to_fit();
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */
	return *this;
}

FMDIndex& FMDIndex::sampleSA(const int32_t* SA, int64_t N, int saSampleRate) {
	assert(N <= INT32_MAX);
	BitStr32 bstr(N);
	SAsampled.clear();
	SAsampled.resize((N + 1) / saSampleRate /* ceil(N / saSampleRate) */ + numGaps()); // enough to hold both peroid sampling and gaps
	for(size_t i = 0, j = 0; i < N; ++i) {
		if(i % saSampleRate == 0 || accessBWT(i) == 0) { /* sample at all null characters */
			bstr.set(i);
			SAsampled[j++] = SA[i];
		}
	}
	SAsampled.shrink_to_fit();
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */
	return *this;
}

FMDIndex& FMDIndex::sampleSA(const int64_t* SA, const DNAseq& bwt, int saSampleRate) {
	const int64_t N = bwt.length();
	BitStr32 bstr(N);
	SAsampled.clear();
	SAsampled.resize((N + 1) / saSampleRate /* ceil(N / saSampleRate) */ + numGaps()); // enough to hold both peroid sampling and gaps
	for(size_t i = 0, j = 0; i < N; ++i) {
		if(i % saSampleRate == 0 || bwt[i] == 0) { /* sample at all null characters */
			bstr.set(i);
			SAsampled[j++] = SA[i];
		}
	}
	SAsampled.shrink_to_fit();
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */
	return *this;
}

FMDIndex& FMDIndex::sampleSA(const int32_t* SA, const DNAseq& bwt, int saSampleRate) {
	const int64_t N = bwt.length();
	assert(N <= INT32_MAX);
	BitStr32 bstr(N);
	SAsampled.clear();
	SAsampled.resize((N + 1) / saSampleRate /* ceil(N / saSampleRate) */ + numGaps()); // enough to hold both peroid sampling and gaps
	for(size_t i = 0, j = 0; i < N; ++i) {
		if(i % saSampleRate == 0 || bwt[i] == 0) { /* sample at all null characters */
			bstr.set(i);
			SAsampled[j++] = SA[i];
		}
	}
	SAsampled.shrink_to_fit();
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */
	return *this;
}

FMDIndex& FMDIndex::buildSA(int saSampleRate) {
	const int64_t N = length();
	DNAseq bwt(N, 0); // temporary storage of 1st-pass accessing values
	{
		/* build the bitstr in the 1st pass */
		BitStr32 bstr(N);
		for(size_t i = 0; i < N; ++i) {
			nt16_t b = accessBWT(i);
			if(i % saSampleRate == 0 || b == 0) /* sample at all null characters */
				bstr.set(i);
			bwt[i] = b;
		}
		SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */
	} /* intermediate bstr will be destoyed */
	SAsampled.resize(SAidx.numOnes()); // enough to hold both peroid sampling and gaps
	/* build SAsampled in the 2nd pass */
#pragma omp parallel for schedule(dynamic)
	for(int64_t i = 0; i < numGaps(); ++i) { // the i-th BWT segment
		int64_t sa = getGapSA(i); // end SA value
		nt16_t b = bwt[i];
		SAsampled[SAidx.rank1(i) - 1] = sa--; // end is always null and sampled
		for(int64_t j = i; b != 0; --sa) {
			j = LF(b, j) - 1; // LF-mapping
			b = bwt[j]; // update b
			if(j % saSampleRate == 0 || b == 0)
				SAsampled[SAidx.rank1(j) - 1] = sa;
		}
	}
	return *this;
}

size_t FMDIndex::getBytes() const {
	return sizeof(B) + sizeof(C) +
			bwtRRR.getBytes() + sizeof(GAParr_t::value_type) * gapSA.size() +
			SAidx.getBytes() + sizeof(SAarr_t::value_type) * SAsampled.size() +
			sizeof(this);
}

int64_t FMDIndex::accessSA(int64_t i) const {
	int64_t dist = 0;
	while(!SAidx.access(i)) {
		i = LF(i) - 1; // backward LF-mapping
		dist++;
	}
	return SAsampled[SAidx.rank1(i) - 1] + dist;
}

vector<GLoc> FMDIndex::locateAllFwd(const DNAseq& pattern) const {
	const size_t L = pattern.length();
	vector<GLoc> locs;
	if(L == 0)
		return locs;

	nt16_t b = pattern.front();
	if(DNAalphabet::isAmbiguous(b))
		return locs;
	int64_t p = C[b];
	int64_t q = C[DNAalphabet::complement(b)];
	int64_t s = C[b + 1] - C[b];

	/* forward search */
    for(int64_t i = 1; i < L && fwdExt(p, q, s, pattern[i]) > 0; ++i)
    	continue;

    for(int64_t j = p; j < p + s; ++j) { // locate fwd locs
    	int64_t start = accessSA(j);
    	locs.push_back(GLoc(start, start + L, -1, GLoc::FWD));
    }
    return locs;
}

vector<GLoc> FMDIndex::locateAllRev(const DNAseq& pattern) const {
	const size_t L = pattern.length();
	vector<GLoc> locs;
	if(L == 0)
		return locs;

	nt16_t b = pattern.back();
	if(DNAalphabet::isAmbiguous(b))
		return locs;
	int64_t p = C[b];
	int64_t q = C[DNAalphabet::complement(b)];
	int64_t s = C[b + 1] - C[b];

	/* backward search */
    for(int64_t i = L - 1; i > 0 && backExt(p, q, s, pattern[i - 1]) > 0; --i)
    	continue;

    for(int64_t j = q; j < q + s; ++j) { // locate rev locs
    	int64_t start = accessSA(j);
    	locs.push_back(GLoc(start, start + pattern.length(), -1, GLoc::REV));
    }
    return locs;
}

int64_t FMDIndex::backExt(int64_t& p, int64_t& q, int64_t& s, nt16_t b) const {
	if(DNAalphabet::isAmbiguous(b))
		return (s = 0);

	array<int64_t, 5> sB;
	array<int64_t, 5> qB;
	int c = DNAalphabet::toInt(b);
	/* calculate new p and s */
	int64_t O = bwtRRR.rank(b, p - 1);
	sB[c] = bwtRRR.rank(b, p + s - 1) - O;

	if(sB[c] <= 0)
		return (s = sB[c]);

	/* update q */
	if(sB[c] != s) {
		sB[0] = bwtRRR.rank(0, p + s - 1) - bwtRRR.rank(0, p - 1);
		for(int j = c + 1; j <= 4 /* T */; ++j) { // search from b + 1 to T
			sB[j] = bwtRRR.rank(DNAalphabet::encode(j), p + s - 1) - bwtRRR.rank(DNAalphabet::encode(j), p - 1);
		}
		/* new range of [q', q' + s' - 1] is a subrange of original [q, q + s] */
		/* devide q + q + s */
		qB[0] = q;
		qB[4] = qB[0] + sB[0];
		for(int j = 4; j > c; --j) // only need to search till b (exclusive)
			qB[j - 1] = qB[j] + sB[j];
		q = qB[c];
	}

	/* update p and s */
	p = C[b] + O;
	s = sB[c];
	return s;
}

BitStr32 FMDIndex::buildInterleavingBS(const FMDIndex& lhs, const FMDIndex& rhs) {
	const int64_t N = lhs.length() + rhs.length();
	BitStr32 bstrM(N);
#pragma omp parallel for schedule(dynamic)
	for(int64_t i = 0; i < lhs.numGaps(); ++i) { // the i-th BWT segment of lhs
		nt16_t b = lhs.accessBWT(i);
		assert(b != 0);
		int64_t RA = rhs.numGaps(); // position on rhs.BWT
#pragma omp critical(WRITE_bstrM)
		bstrM.set(i + RA);
		for(int64_t j = i; b != 0; b = lhs.accessBWT(j)) { // j is position on lhs.BWT
			j = lhs.LF(b, j) - 1; // LF-mapping on lhs
			RA = rhs.LF(b, RA - 1); // LF-mapping on rhs
#pragma omp critical(WRITE_bstrM)
			bstrM.set(j + RA);
		}
	}
	return bstrM;
}

DNAseq FMDIndex::mergeBWT(const FMDIndex& lhs, const FMDIndex& rhs, const BitStr32& bstrM) {
	const int64_t N = lhs.length() + rhs.length();
	assert(N == bstrM.length());
	DNAseq bwtM(N, 0); // merged BWT
	for(size_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM[k] = bstrM.test(k) ? lhs.accessBWT(i++) : rhs.accessBWT(j++);
	return bwtM;
}

DNAseq FMDIndex::mergeBWT(const FMDIndex& lhs, const FMDIndex& rhs, const BitSeq& bsM) {
	const int64_t N = lhs.length() + rhs.length();
	assert(N == bsM.length());
	DNAseq bwtM(N, 0); // merged BWT
#pragma omp parallel for
	for(size_t k = 0; k < N; ++k)
		bwtM[k] = bsM.test(k) ? lhs.accessBWT(bsM.rank1(k) - 1) : rhs.accessBWT(bsM.rank0(k) - 1);
	return bwtM;
}

FMDIndex::GAParr_t FMDIndex::mergeGap(const FMDIndex& lhs, const FMDIndex& rhs) {
	const int64_t N1 = lhs.numGaps();
	const int64_t N2 = rhs.numGaps();
	const int64_t N = N1 + N2;
	GAParr_t gapM(N, 0); // init with 0s
	/* copy rhs gaps with shift */
	std::transform(rhs.gapSA.begin(), rhs.gapSA.end(), gapM.begin(),
			[&] (GAParr_t::value_type pos) { return pos + lhs.length(); } );
	/* copy lhs gaps */
	std::copy(lhs.gapSA.begin(), lhs.gapSA.end(), gapM.begin() + N2);
	return gapM;
}

FMDIndex& FMDIndex::mergeCount(const FMDIndex& other) {
	for(int64_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
		B[i] += other.B[i];
		C[i] += other.C[i];
	}
	return *this;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
