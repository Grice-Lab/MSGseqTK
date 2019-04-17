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
#include "FMDIndex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace EGriceLab {
namespace MSGseqTK {
using std::vector;
using EGriceLab::libSDS::BitStr32;

FMDIndex& FMDIndex::build(const DNAseq& seq) {
	assert(seq.back() == DNAalphabet::GAP_BASE);
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
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
	bwt.save(out);
	out.write((const char*) &keepSA, sizeof(bool));
	if(keepSA) {
		size_t nSAsampled = SAsampled.size();
		assert(nSAsampled == SAbit.numOnes());
		out.write((const char*) &nSAsampled, sizeof(size_t));
		out.write((const char*) SAsampled.data(), sizeof(int64_t) * nSAsampled);
		SAbit.save(out);
	}
	return out;
}

istream& FMDIndex::load(istream& in) {
	in.read((char*) B.data(), B.size() * sizeof(BCarray_t::value_type));
	in.read((char*) C.data(), C.size() * sizeof(BCarray_t::value_type));
	bwt.load(in);
	in.read((char*) &keepSA, sizeof(bool));
	if(keepSA) {
		size_t nSAsampled;
		in.read((char*) &nSAsampled, sizeof(size_t));
		SAsampled.resize(nSAsampled);
		in.read((char*) SAsampled.data(), sizeof(int64_t) * nSAsampled);
		SAbit.load(in);
		assert(nSAsampled == SAbit.numOnes());
	}
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

	BitStr32 bstrM = buildInterleavingBS(*this, other); // build interleaving bitvector
	DNAseq bwtM = mergeBWT(this->bwt, other.bwt, bstrM); // merge BWT
	bstrM.reset(); // clear RAM
    bwt = WaveletTreeRRR(bwtM, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE); // update bwtRRR

	/* merging counts after bwt updated */
	for(int64_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
		B[i] += other.B[i];
		C[i] += other.C[i];
	}

    /* reset SA */
	SAbit.reset();

    if(keepSA)
    	buildSA(bwtM);

	return *this;
}

DNAseq FMDIndex::getBWT() const {
	DNAseq bwtSeq;
	bwtSeq.reserve(length());
	for(int64_t i = 0; i < length(); ++i)
		bwtSeq.push_back(bwt.access(i));
	return bwtSeq;
}

DNAseq FMDIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	DNAseq seq;
	seq.reserve(length());
	for(int64_t i = 0; i < C[0 + 1]; ++i) { // i-th pass of LF-mapping
		seq.push_back(0);
		uint8_t b;
		for(int64_t j = i; (b = bwt.access(j)) != 0; j = LF(b, j) - 1 /* LF-mapping */)
			seq.push_back(b);
	}
	assert(seq.length() == length());
	std::reverse(seq.begin(), seq.end()); // LF-transverse is on reverse (suffix) order
	return seq;
}

void FMDIndex::buildBWT(const DNAseq& seq) {
	assert(seq.back() == DNAalphabet::GAP_BASE);
	const size_t N = seq.length();
	/* construct SA */
    int64_t* SA = new int64_t[N];
    int64_t errn = divsufsort((const uint8_t*) (seq.c_str()), (saidx_t*) SA, N);
	if(errn != 0)
		throw std::runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build bwt and sample bitstr */
	DNAseq bwtSeq;
	bwtSeq.reserve(N);
	for(int64_t* sa = SA; sa < SA + N; ++sa)
		bwtSeq.push_back(*sa == 0 ? 0 : seq[*sa - 1]);

	/* construct BWTRRR */
    bwt = WaveletTreeRRR(bwtSeq, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);

	if(keepSA)
		buildSA(SA, bwtSeq);
	delete[] SA;
}

void FMDIndex::buildSA() {
	assert(isInitiated());
	const int64_t N = length();
	/* build a BitStr in the 1st pass */
	BitStr32 bstr(N);
#pragma omp parallel for
	for(int64_t i = 0; i < N; ++i) {
		if(i % SA_SAMPLE_RATE == 0 || bwt.access(i) == 0) /* sample at all null characters */
			bstr.set(i);
	}
	SAbit = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(SAbit.numOnes()); /* sample at on bits */
	int64_t k = N;
	for(int64_t i = 0; i < C[0 + 1]; ++i) { // i-th pass of LF-mapping loop
		int64_t j = i;
		uint8_t b;
		do {
			b = bwt.access(j);
			if(bstr.test(j))
				SAsampled[SAbit.rank1(j) - 1] = k - 1;
			j = LF(b, j) - 1; // LF-mapping
			k--;
		}
		while(b != 0);
	}
	assert(k == 0);
}

void FMDIndex::buildSA(const DNAseq& bwtSeq) {
	const int64_t N = length();
	assert(bwtSeq.length() == N);
	/* build the bitstr by sampling bwtSeq */
	BitStr32 bstr(N);
	for(size_t i = 0; i < N; ++i)
		if(i % SA_SAMPLE_RATE == 0 || bwtSeq[i] == 0) /* sample at all null characters */
			bstr.set(i);

	SAbit = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(SAbit.numOnes()); /* sample at on bits */
	int64_t k = N; // position on seq, start from 1 after end
	for(int64_t i = 0; i < C[0 + 1]; ++i) { // the i-th null segment
		int64_t j = i;
		uint8_t b;
		do {
			b = bwtSeq[j];
			if(bstr.test(j))
				SAsampled[SAbit.rank1(j) - 1] = k - 1;
			j = LF(b, j) - 1; // LF-mapping
			k--;
		}
		while(b != 0);
	}
	assert(k == 0);
}

void FMDIndex::buildSA(const int64_t* SA, const DNAseq& bwtSeq) {
	const int64_t N = length();
	assert(bwtSeq.length() == N);
	/* build the bitstr by sampling SA direction */
	BitStr32 bstr(N);
	for(size_t i = 0; i < N; ++i)
		if(i % SA_SAMPLE_RATE == 0 || bwtSeq[i] == 0) /* sample at all null characters */
			bstr.set(i);

	SAbit = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(SAbit.numOnes()); /* sample at on bits */
#pragma omp parallel for
	for(int64_t i = 0; i < N; ++i) {
		if(bstr.test(i)) {
			SAsampled[SAbit.rank1(i) - 1] = SA[i]; // no lock needed since independent access guaranteed
		}
	}
}

int64_t FMDIndex::accessSA(int64_t i) const {
	int64_t dist = 0;
	while(!SAbit.access(i)) {
		i = LF(i) - 1; // backward LF-mapping
		dist++;
	}
	return SAsampled[SAbit.rank1(i) - 1] + dist;
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

FMDIndex::FMDIndex(const BCarray_t& B, const BCarray_t& C, const DNAseq& bwtSeq, bool keepSA)
: B(B), C(C), bwt(WaveletTreeRRR(bwtSeq, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE)), keepSA(keepSA)
{
	/* build SA using customized bitstr */
	if(keepSA)
		buildSA(bwtSeq);
}

FMDIndex operator+(const FMDIndex& lhs, const FMDIndex& rhs) {
	if(!rhs.isInitiated())
		return lhs;
	if(!lhs.isInitiated())
		return rhs;

	const int64_t N1 = lhs.length();
	const int64_t N2 = rhs.length();
	const int64_t N = N1 + N2;
	FMDIndex::BCarray_t BMerged;
	FMDIndex::BCarray_t CMerged;

	BitStr32 bstrM = FMDIndex::buildInterleavingBS(lhs, rhs); // build interleaving bitvector
	DNAseq bwtM = FMDIndex::mergeBWT(lhs.bwt, rhs.bwt, bstrM); // merge BWT
	bstrM.reset(); // clear RAM

	/* merge counts */
	for(int64_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
		BMerged[i] = lhs.B[i] + rhs.B[i];
		CMerged[i] = lhs.C[i] + rhs.C[i];
	}

	return FMDIndex(BMerged, CMerged, bwtM, lhs.keepSA || rhs.keepSA /* keep SA if any of the two operands keepSA */);
}

int64_t FMDIndex::backExt(int64_t& p, int64_t& q, int64_t& s, uint8_t b) const {
	if(!DNAalphabet::isBasic(b))
		return 0;
	int64_t pN; // fwd strand backExt
	BCarray_t qB, sB;
	qB.fill(0);
	sB.fill(0);
	/* calculate new p and s */
	int64_t O = bwt.rank(b, p - 1);
	pN = C[b] + O;
	sB[b] = bwt.rank(b, p + s - 1) - O;

	/* update q if s changes */
	if(sB[b] != s) {
		sB[0] = bwt.rank(0, p + s - 1) - bwt.rank(0, p - 1);
		for(nt16_t i = b + 1; i <= DNAalphabet::NT16_MAX; ++i) { // search from b + 1
			if(DNAalphabet::isBasic(i))
				sB[i] = bwt.rank(i, p + s - 1) - bwt.rank(i, p - 1);
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

BitStr32 FMDIndex::buildInterleavingBS(const FMDIndex& fmdidx1, const FMDIndex& fmdidx2) {
	const size_t N = fmdidx1.length() + fmdidx2.length();
	BitStr32 bstrM(N);
#pragma omp parallel for
	for(int64_t i = 0; i < fmdidx1.C[0 + 1]; ++i) { // i-th pass LF-mapping on lhs
		int64_t j = i;
		int64_t RA = fmdidx2.C[0 + 1];
		uint8_t b;
		do {
#pragma omp critical(WRITE_bstrM)
			bstrM.set(j + RA);
			b = fmdidx1.bwt.access(j);
			/* LF-mapping */
			RA = fmdidx2.LF(b, RA - 1);
			j = fmdidx1.LF(b, j) - 1;
		}
		while(b != 0);
	}
	return bstrM;
}

DNAseq FMDIndex::mergeBWT(const WaveletTreeRRR& bwt1, const WaveletTreeRRR& bwt2, const BitStr32& bstr) {
	const size_t N = bstr.length();
	assert(N == bwt1.length() + bwt2.length());
	DNAseq bwtM(N, 0);
#ifndef _OPENMP
	for(size_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM[k] = bstr.test(k) ? bwt1.access(i++) : bwt2.access(j++);
#else
	int nThreads = omp_get_num_threads();
	if(nThreads == 1) { // no parallelzation needed
		for(size_t i = 0, j = 0, k = 0; k < N; ++k)
			bwtM[k] = bstr.test(k) ? bwt1.access(i++) : bwt2.access(j++);
	}
	else {
		BitSeqRRR bs(bstr); // aux BitSeq for rank access
#pragma omp parallel for
		for(size_t k = 0; k < N; ++k)
			bwtM[k] = bstr.test(k) ? bwt1.access(bs.rank1(k) - 1) : bwt2.access(bs.rank0(k) - 1);
	}
#endif
	return bwtM;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
