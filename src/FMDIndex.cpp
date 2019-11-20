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

FMDIndex::FMDIndex(const DNAseq& seq, bool keepSA, int saSampleRate) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");
	buildCounts(seq); // build count
	const int64_t* SA = buildBWT(seq); // build BWT
	buildGap(SA); // build GAP
	if(keepSA)
		buildSA(SA, saSampleRate);
	else
		clearSA();
	delete[] SA; // delete temporary
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
	StringUtils::saveString(gapSA, out);
	SAidx.save(out);
	StringUtils::saveString(SAsampled, out);
	return out;
}

istream& FMDIndex::load(istream& in) {
	in.read((char*) B.data(), B.size() * sizeof(BCarray_t::value_type));
	in.read((char*) C.data(), C.size() * sizeof(BCarray_t::value_type));
	StringUtils::loadString(bwt, in);
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

FMDIndex& FMDIndex::operator+=(const FMDIndex& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

    /* merge gap info */
    gapSA = mergeGap(*this, other);

	/* merge BWTs */
	bwt = mergeBWT(*this, other, buildInterleavingBS(*this, other)); // update uncompressed BWT
    bwtRRR = WaveletTreeRRR(bwt, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR

	/* merge counts */
	for(int64_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
		B[i] += other.B[i];
		C[i] += other.C[i];
	}

	/* reset SA but do not update */
	clearSA();
	return *this;
}

DNAseq FMDIndex::getBWT() const {
	if(hasBWT())
		return bwt;
	DNAseq bwtSeq;
	bwtSeq.reserve(length());
	for(int64_t i = 0; i < length(); ++i)
		bwtSeq.push_back(bwtRRR.access(i));
	return bwtSeq;
}

DNAseq FMDIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	const int64_t N = length();
	DNAseq seq(N, 0);
#pragma omp parallel for schedule(dynamic, 4)
	for(int64_t i = 0; i < numGaps(); ++i) { // i-th BWT segment
		int64_t sa = getGapSA(i); // end SA value
		seq[sa] = 0;
		nt16_t b = accessBWT(i);
		sa--;
		for(int64_t j = i; b != 0; b = accessBWT(j)) {
			seq[sa] = b;
			j = LF(b, j) - 1; // LF-mapping
			sa--;
		}
	}
//	assert(seq.length() == N);
	return seq;
}

int64_t* FMDIndex::buildBWT(const DNAseq& seq) {
	assert(seq.back() == 0);
	const size_t N = seq.length();
	/* construct SA */
	int64_t* SA = new int64_t[N];
    int64_t errn = divsufsort((const uint8_t*) seq.c_str(), (saidx_t*) SA, N); // string.c_str guarantee a null terminal
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

FMDIndex& FMDIndex::buildSA(int saSampleRate) {
	const int64_t N = length();
	assert(bwt.length() == N);
	/* build the bitstr in the 1st pass */
	BitStr32 bstr(N);
	for(size_t i = 0; i < N; ++i) {
		if(i % saSampleRate == 0 || bwt[i] == 0) /* sample at all null characters */
			bstr.set(i);
	}
	SAidx = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */

	SAsampled.clear();
	SAsampled.resize(SAidx.numOnes()); // enough to hold both peroid sampling and gaps
	/* build SAsampled in the 2nd pass */
#pragma omp parallel for schedule(dynamic, 4)
	for(int64_t i = 0; i < numGaps(); ++i) { // the i-th BWT segment
		int64_t sa = getGapSA(i); // end SA value
		nt16_t b = bwt[i];
		SAsampled[SAidx.rank1(i) - 1] = sa--; // end is always null and sampled
		for(int64_t j = i; b != 0; b = bwt[j]) {
			j = LF(b, j) - 1; // LF-mapping
			if(SAidx.test(j))
				SAsampled[SAidx.rank1(j) - 1] = sa;
			sa--;
		}
	}
	return *this;
}

void FMDIndex::buildGap(const int64_t* SA) {
	/* build gapSA */
	gapSA.resize(numGaps());
	for(int64_t i = 0; i < numGaps(); ++i)
		gapSA[i] = SA[i];
}

void FMDIndex::buildSA(const int64_t* SA, int saSampleRate) {
	const int64_t N = length();
	assert(bwt.length() == N);
	/* build the bitstr by sampling SA direction */
	BitStr32 bstr(N);
	SAsampled.clear();
	SAsampled.reserve(N / saSampleRate + numGaps()); // enough to hold both peroid sampling and gaps
	for(size_t i = 0; i < N; ++i) {
		if(i % saSampleRate == 0 || bwt[i] == 0) { /* sample at all null characters */
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
	if(DNAalphabet::isAmbiguous(b)) {
		s = 0;
		return s;
	}
	int64_t pN; // fwd strand backExt
	BCarray_t qB {};
	BCarray_t sB {};
//	qB.fill(0);
//	sB.fill(0);
	/* calculate new p and s */
	int64_t O = bwtRRR.rank(b, p - 1);
	pN = C[b] + O;
	sB[b] = bwtRRR.rank(b, p + s - 1) - O;

	/* update q */
	if(sB[b] != s) {
		sB[0] = bwtRRR.rank(0, p + s - 1) - bwtRRR.rank(0, p - 1);
		for(nt16_t i : { DNAalphabet::A, DNAalphabet::C, DNAalphabet::G, DNAalphabet::T }) { // search from b + 1
			if(i > b)
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

FMDIndex::GAParr_t FMDIndex::mergeGap(const FMDIndex& lhs, const FMDIndex& rhs) {
	const int64_t N1 = lhs.numGaps();
	const int64_t N2 = rhs.numGaps();
	GAParr_t gapM(N1 + N2, 0); // init with 0s
	/* copy rhs gaps with shift */
	std::transform(rhs.gapSA.begin(), rhs.gapSA.end(), gapM.begin(),
			[&] (GAParr_t::value_type pos) { return pos + lhs.length(); } );
	/* copy lhs gaps */
	std::copy(lhs.gapSA.begin(), lhs.gapSA.end(), gapM.begin() + N2);
	return gapM;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
