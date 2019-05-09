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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace EGriceLab {
namespace MSGseqTK {
using std::vector;
using EGriceLab::libSDS::BitStr32;
using EGriceLab::libSDS::BitSeqGGMN; // useful for temporary and fast BitSeq solution

FMDIndex::FMDIndex(const DNAseq& seq, bool keepSA) {
	assert(seq.back() == DNAalphabet::GAP_BASE);
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");
	const size_t N = seq.length();
	buildCounts(seq); // build count
	const int64_t* SA = buildBWT(seq); // build BWT
	if(keepSA)
		buildSA(SA);
	else
		clearSA();
	delete[] SA; // delete temporary
}

void FMDIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::value_type b : seq)
		B[b]++;

	/* calculate cumulative counts */
    int64_t S = 0;
    for(nt16_t i = 0; i < DNAalphabet::SIZE; ++i) {
    	C[i] = S;
    	S += B[i];
    }
    if(!isBiDirectional()) {
    	std::cerr << "input seq is not bi-directional" << std::endl;
    	abort();
    }
    if(getExtBaseCount() > 0) {
    	std::cerr << "input seq does not allow IUPAC-extended bases" << std::endl;
    	abort();
    }
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
	StringUtils::saveString(bwt, out);
	bwtRRR.save(out);
	SAidx.save(out);
	StringUtils::saveString(SAsampled, out);
//	StringUtils::saveString(SAgap, out);
	return out;
}

istream& FMDIndex::load(istream& in) {
	in.read((char*) B.data(), B.size() * sizeof(BCarray_t::value_type));
	in.read((char*) C.data(), C.size() * sizeof(BCarray_t::value_type));
	StringUtils::loadString(bwt, in);
	bwtRRR.load(in);
	SAidx.load(in);
	StringUtils::loadString(SAsampled, in);
//	StringUtils::loadString(SAgap, in);
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

FMDIndex& FMDIndex::operator+=(const FMDIndex& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

	bwt = mergeBWT(*this, other, buildInterleavingBS(*this, other)); // update uncompressed BWT
    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR

	/* merging counts after bwt updated */
	for(int64_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
		B[i] += other.B[i];
		C[i] += other.C[i];
	}

	/* reset SA but do not update */
	clearSA();
	return *this;
}

DNAseq FMDIndex::getBWT() const {
	DNAseq bwtSeq;
	bwtSeq.reserve(length());
	for(int64_t i = 0; i < length(); ++i)
		bwtSeq.push_back(bwtRRR.access(i));
	return bwtSeq;
}

DNAseq FMDIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	DNAseq seq;
	seq.reserve(length());
	for(int64_t i = 0; i < B[0]; ++i) { // i-th pass of LF-mapping
		seq.push_back(0);
		nt16_t b;
		for(int64_t j = i; (b = accessBWT(j)) != 0; j = LF(b, j) - 1 /* LF-mapping */)
			seq.push_back(b);
	}
	assert(seq.length() == length());
	std::reverse(seq.begin(), seq.end()); // LF-transverse is on reverse (suffix) order
	return seq;
}

int64_t* FMDIndex::buildBWT(const DNAseq& seq) {
	assert(seq.back() == DNAalphabet::GAP_BASE);
	const size_t N = seq.length();
	/* construct SA */
	int64_t* SA = new int64_t[N];
    int64_t errn = divsufsort((const uint8_t*) seq.data(), (saidx_t*) SA, N);
	if(errn != 0)
		throw std::runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build uncompressed bwt */
	bwt.clear();
	bwt.reserve(N);
	for(const int64_t* sa = SA; sa < SA + N; ++sa)
		bwt.push_back(*sa == 0 ? 0 : seq[*sa - 1]);

	/* build BWTRRR */
    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);
    return SA;
}

FMDIndex& FMDIndex::buildSA() {
	const int64_t N = length();
	assert(bwt.length() == N);
//	int64_t* SA = new int64_t[N];
//	int64_t k = N; // position on seq/SA, start from 1 pass SA end
//	for(int64_t i = 0; i < B[0]; ++i) { // the i-th null segment
//		nt16_t b = bwt[i];
//		SA[i] = --k;
//		for(int64_t j = i; b != 0; b = bwt[j]) {
////			std::cerr << "i: " << i << " j: " << j << " k: " << k << std::endl;
//			j = LF(b, j) - 1; // LF-mapping
//			SA[j] = --k;
//		}
//	}
//	assert(k == 0);
//	buildSA(SA);
//	delete[] SA;

	/* build the bitstr by sampling bwtSeq */
	BitStr32 bstr(N);
	for(size_t i = 0; i < N; ++i) {
		if(i % SA_SAMPLE_RATE == 0 || bwt[i] == 0) /* sample at all null characters */
			bstr.set(i);
	}
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* update SAidx */

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(SAidx.numOnes()); /* sample at on bits */
	int64_t k = N; // position on seq/SA, start from 1 pass SA end
	for(int64_t i = 0; i < B[0]; ++i) { // the i-th null segment
		nt16_t b = bwt[i];
		SAsampled[SAidx.rank1(i) - 1] = --k; // always sample at null
		for(int64_t j = i; b != 0; b = bwt[j]) {
			//	std::cerr << "i: " << i << " j: " << j << " k: " << k << " b: " << DNAalphabet::decode(b) << std::endl;
			j = LF(b, j) - 1; // LF-mapping
			if(bstr.test(j))
				SAsampled[SAidx.rank1(j) - 1] = k - 1;
			k--;
		}
	}
	assert(k == 0);
	return *this;
}

FMDIndex& FMDIndex::buildSA(const vector<size_t>& gapLoc) {
	assert(gapLoc.size() == B[0]);
	assert(gapLoc.back() == length() - 1); // last base must be a GAP_BASE
	const int64_t N = length();
	assert(bwt.length() == N);
	/* build the bitstr by sampling bwtSeq */
	BitStr32 bstr(N);
	for(size_t i = 0; i < N; ++i) {
		if(i % SA_SAMPLE_RATE == 0 || bwt[i] == 0) /* sample at all null characters */
			bstr.set(i);
	}
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* update SAidx */

	/* parallel build SAsampled in the 2nd pass */
	SAsampled.resize(SAidx.numOnes()); /* sample at on bits */
#pragma omp parallel for schedule(dynamic, 4)
	for(int64_t i = 0; i < B[0]; ++i) { // the i-th null segment
		int64_t j = i; // position on BWT
		int64_t k = gapLoc[B[0] - i - 1]; // search is backward
		nt16_t b = bwt[i];
		SAsampled[SAidx.rank1(i) - 1] = k--;
		for(int64_t j = i; b != 0; b = bwt[j]) {
			j = LF(b, j) - 1; // LF-mapping
			if(bstr.test(j))
				SAsampled[SAidx.rank1(j) - 1] = k--;
		}
	}
	return *this;
}

void FMDIndex::buildSA(const int64_t* SA) {
	const int64_t N = length();
	assert(bwt.length() == N);
	/* build the bitstr by sampling SA direction */
	BitStr32 bstr(N);
	SAsampled.clear();
	SAsampled.reserve(N / SA_SAMPLE_RATE + B[0]); // enough to hold both peroid sampling and gaps
	for(size_t i = 0; i < N; ++i) {
		if(i % SA_SAMPLE_RATE == 0 || bwt[i] == 0) { /* sample at all null characters */
			bstr.set(i);
			SAsampled.push_back(SA[i]);
		}
	}
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */
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

	size_t i = 0;
	nt16_t b = pattern[i];
	int64_t p = C[b];
	int64_t q = C[DNAalphabet::complement(b)];
	int64_t s = C[b + 1] - C[b];

	/* forward search */
    for(i = 1; i < L && fwdExt(p, q, s, pattern[i]) > 0; ++i)
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

	size_t i = L;
	nt16_t b = pattern[i - 1];
	int64_t p = C[b];
	int64_t q = C[DNAalphabet::complement(b)];
	int64_t s = C[b + 1] - C[b];

	/* backward search */
    for(i = L - 1; i > 0 && backExt(p, q, s, pattern[i - 1]) > 0; --i)
    	continue;

    for(int64_t j = q; j < q + s; ++j) { // locate rev locs
    	int64_t start = accessSA(j);
    	locs.push_back(GLoc(start, start + pattern.length(), -1, GLoc::REV));
    }
    return locs;
}

int64_t FMDIndex::backExt(int64_t& p, int64_t& q, int64_t& s, nt16_t b) const {
	if(!DNAalphabet::isBasic(b)) {
		s = 0;
		return s;
	}
	int64_t pN; // fwd strand backExt
	BCarray_t qB, sB;
	qB.fill(0);
	sB.fill(0);
	/* calculate new p and s */
	int64_t O = bwtRRR.rank(b, p - 1);
	pN = C[b] + O;
	sB[b] = bwtRRR.rank(b, p + s - 1) - O;

	/* update q if s changes */
	if(sB[b] != s) {
		sB[0] = bwtRRR.rank(0, p + s - 1) - bwtRRR.rank(0, p - 1);
		for(nt16_t i = b + 1; i <= DNAalphabet::NT16_MAX; ++i) { // search from b + 1
			if(DNAalphabet::isBasic(i))
				sB[i] = bwtRRR.rank(i, p + s - 1) - bwtRRR.rank(i, p - 1);
		}
		/* new range of [q', q' + s' - 1] is a subrange of original [q, q + s] */
		/* devide q + q + s */
		qB[0] = q;
		qB[DNAalphabet::T] = qB[0] + sB[0];
		for(nt16_t i = DNAalphabet::T; i > b; --i) // only need to search till b + 1
			qB[i - 1] = qB[i] + sB[i];
		q = qB[b];
	}

	/* update p and s */
	p = pN;
	s = sB[b];
	return s;
}

BitStr32 FMDIndex::buildInterleavingBS(const FMDIndex& lhs, const FMDIndex& rhs) {
	const size_t N = lhs.length() + rhs.length();
	BitStr32 bstrM(N);
#pragma omp parallel for schedule(dynamic, 4)
	for(int64_t i = 0; i < lhs.B[0]; ++i) { // i-th pass LF-mapping on lhs
		int64_t j = i; // position on lhs.BWT
		nt16_t b = lhs.accessBWT(i);
		int64_t RA = rhs.B[0]; // position on rhs.BWT
#pragma omp critical(WRITE_bstrM)
		bstrM.set(i + RA);
		for(int64_t j = i; b != 0; b = lhs.accessBWT(j)) {
			j = lhs.LF(b, j) - 1; // LF-mapping
			RA = rhs.LF(b, RA - 1);
#pragma omp critical(WRITE_bstrM)
			bstrM.set(j + RA);
		}
	}
	return bstrM;
}

DNAseq FMDIndex::mergeBWT(const FMDIndex& lhs, const FMDIndex& rhs, const BitStr32& bstrM) {
	const size_t N = lhs.length() + rhs.length();
	assert(N == bstrM.length());
	DNAseq bwtM(N, 0); // merged BWT
#ifndef _OPENMP
	for(size_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM[k] = bstrM.test(k) ? lhs.accessBWT(i++) : rhs.accessBWT(j++);
#else
	int nThreads = 1;
#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}
	if(1 == nThreads) { // no parallelzation needed
		for(size_t i = 0, j = 0, k = 0; k < N; ++k)
			bwtM[k] = bstrM.test(k) ? lhs.accessBWT(i++) : rhs.accessBWT(j++);
	}
	else {
		BitSeqGGMN bsM(bstrM);
#pragma omp parallel for
		for(size_t k = 0; k < N; ++k) {
			bwtM[k] = bstrM.test(k) ? lhs.accessBWT(bsM.rank1(k) - 1) : rhs.accessBWT(bsM.rank0(k) - 1);
		}
	}
#endif
	return bwtM;
}

DNAseq FMDIndex::mergeBWT(const FMDIndex& lhs, const FMDIndex& rhs, BitStr32&& bstrM) {
	const size_t N = lhs.length() + rhs.length();
	assert(N == bstrM.length());
	DNAseq bwtM(N, 0); // merged BWT
#ifndef _OPENMP
	for(size_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM[k] = bstrM.test(k) ? lhs.accessBWT(i++) : rhs.accessBWT(j++);
#else
	int nThreads = 1;
#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}
	if(1 == nThreads) { // no parallelzation needed
		for(size_t i = 0, j = 0, k = 0; k < N; ++k)
			bwtM[k] = bstrM.test(k) ? lhs.accessBWT(i++) : rhs.accessBWT(j++);
	}
	else {
		BitSeqGGMN bsM(bstrM);
#pragma omp parallel for
		for(size_t k = 0; k < N; ++k) {
			bwtM[k] = bstrM.test(k) ? lhs.accessBWT(bsM.rank1(k) - 1) : rhs.accessBWT(bsM.rank0(k) - 1);
		}
	}
#endif
	return bwtM;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
