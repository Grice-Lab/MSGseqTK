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

const saidx_t FMDIndex::totalBases() const {
	saidx_t N = 0;
	for(saidx_t i = 0; i < B.size(); ++i)
		N += B[i];
	return N;
}

FMDIndex& FMDIndex::build(const DNAseq& seq) {
	if(seq.length() > MAX_LENGTH)
		throw std::length_error("DNAseq length exceeding the max allowed length");

	buildCounts(seq);
	buildBWT(seq);

	return *this;
}

void FMDIndex::buildCounts(const DNAseq& seq) {
	for(DNAseq::value_type b : seq)
		B[b]++;
	B[0]++; // count null terminal

	/* calculate cumulative counts */
    saidx_t S = 0;
    for(nt16_t i = 0; i < DNAalphabet::SIZE; ++i) {
    	C[i] = S;
    	S += B[i];
    }
    if(!isBiDirectional()) {
    	std::cerr << "input seq is not bi-directional" << std::endl;
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
		out.write((const char*) SAsampled.data(), sizeof(saidx_t) * nSAsampled);
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
		in.read((char*) SAsampled.data(), sizeof(saidx_t) * nSAsampled);
		SAbit.load(in);
		assert(nSAsampled == SAbit.numOnes());
	}
	return in;
}

saidx_t FMDIndex::count(const DNAseq& pattern) const {
	if(pattern.empty())
		return 0;

	DNAseq::const_reverse_iterator b = pattern.rbegin(); /* could be the null terminal */
	saidx_t p = C[*b];
	saidx_t q = C[DNAalphabet::complement(*b)];
	saidx_t s = C[*b + 1] - C[*b];

	/* backward search */
    while(s > 0 && ++b < pattern.rend()) {
    	backExt(p, q, s, *b);
    }
    return std::max<int64_t>(s, 0);
}

FMDIndex& FMDIndex::operator+=(const FMDIndex& other) {
	if(!other.isInitiated()) /* cannot merge */
		return *this;
	if(!isInitiated()) {
		*this = other;
		return *this;
	}

	/* build interleaving bitvector between *this and other */
	const saidx_t N1 = length();
	const saidx_t N2 = other.length();
	const saidx_t N = N1 + N2;

	/* build RA and interleaving bitvector */
	BitStr32 bstr(N);
#pragma omp parallel for
	for(saidx_t i = 0; i < C[0 + 1]; ++i) { // i-th pass of LF-mapping
		saidx_t j = i;
		saidx_t RA = other.C[0 + 1];
		sauchar_t b;
	  do {
#pragma omp critical(WRITE_bstr)
			bstr.set(j + RA);
			b = bwt.access(j);
			/* LF-mapping */
			RA = other.LF(b, RA - 1);
			j = LF(b, j) - 1;
		}
		while(b != 0);
	}

	/* build merbed BWT */
	DNAseq bwtM;
	bwtM.reserve(N);
	for(saidx_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM.push_back(bstr.test(k) ? bwt.access(i++) : other.bwt.access(j++));
	/* update bwtRRR */
    bwt = WaveletTreeRRR(bwtM, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);

	/* merging B[] and C[] */
	for(saidx_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
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
	for(saidx_t i = 0; i < length(); ++i)
		bwtSeq.push_back(bwt.access(i));
	return bwtSeq;
}

DNAseq FMDIndex::getSeq() const {
	/* get Seq by LF-mapping transverse */
	DNAseq seq;
	seq.reserve(length());
	for(saidx_t i = 0; i < C[0 + 1]; ++i) { // i-th pass of LF-mapping
		seq.push_back(0);
		sauchar_t b;
		for(saidx_t j = i; (b = bwt.access(j)) != 0; j = LF(b, j) - 1 /* LF-mapping */)
			seq.push_back(b);
	}
	assert(seq.length() == length());
	std::reverse(seq.begin(), seq.end()); // LF-transverse is on reverse (suffix) order
	return seq;
}

void FMDIndex::buildBWT(const DNAseq& seq) {
	const size_t N =  seq.length() + 1;
	/* construct SA */
    saidx_t* SA = new saidx_t[N];
    saidx_t errn = divsufsort((const sauchar_t*) (seq.c_str()), SA, N);
	if(errn != 0)
		throw std::runtime_error("Error: Cannot build suffix-array on DNAseq");

	/* build bwt and sample bitstr */
	DNAseq bwtSeq;
	bwtSeq.reserve(N);
	for(saidx_t* sa = SA; sa < SA + N; ++sa)
		bwtSeq.push_back(*sa == 0 ? 0 : seq[*sa - 1]);

	/* construct BWTRRR */
    bwt = WaveletTreeRRR(bwtSeq, 0, DNAalphabet::NT16_MAX, RRR_SAMPLE_RATE);

	if(keepSA)
		buildSA(SA, bwtSeq);
	delete[] SA;
}

void FMDIndex::buildSA() {
	assert(isInitiated());
	const saidx_t N = length();
	/* build a BitStr in the 1st pass */
	BitStr32 bstr(N);
#pragma omp parallel for
	for(saidx_t i = 0; i < N; ++i) {
		if(i % SA_SAMPLE_RATE == 0 || bwt.access(i) == 0) /* sample at all null characters */
			bstr.set(i);
	}
	SAbit = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(SAbit.numOnes()); /* sample at on bits */
	saidx_t k = N;
	for(saidx_t i = 0; i < C[0 + 1]; ++i) { // i-th pass of LF-mapping loop
		saidx_t j = i;
		sauchar_t b;
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
	const saidx_t N = length();
	assert(bwtSeq.length() == N);
	/* build the bitstr by sampling bwtSeq */
	BitStr32 bstr(N);
	for(size_t i = 0; i < N; ++i)
		if(i % SA_SAMPLE_RATE == 0 || bwtSeq[i] == 0) /* sample at all null characters */
			bstr.set(i);

	SAbit = BitSeqRRR(bstr, RRR_SAMPLE_RATE); /* reset the SAbit */

	/* build SAsampled in the 2nd pass */
	SAsampled.resize(SAbit.numOnes()); /* sample at on bits */
	saidx_t k = N; // position on seq, start from 1 after end
	for(saidx_t i = 0; i < C[0 + 1]; ++i) { // the i-th null segment
		saidx_t j = i;
		sauchar_t b;
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

void FMDIndex::buildSA(const saidx_t* SA, const DNAseq& bwtSeq) {
	const saidx_t N = length();
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
	for(saidx_t i = 0; i < N; ++i) {
		if(bstr.test(i)) {
			SAsampled[SAbit.rank1(i) - 1] = SA[i];
		}
	}
}

saidx_t FMDIndex::accessSA(saidx_t i) const {
	saidx_t dist = 0;
	while(!SAbit.access(i)) {
		i =  LF(i) - 1; // backward LF-mapping
		dist++;
	}
	return SAsampled[SAbit.rank1(i) - 1] + dist;
}

vector<GLoc> FMDIndex::locateAllFwd(const DNAseq& pattern) const {
	vector<GLoc> locs;
	if(pattern.empty())
		return locs;

	nt16_t b = pattern.front();
	saidx_t p = C[b];
	saidx_t q = C[DNAalphabet::complement(b)];
	saidx_t s = C[b + 1] - C[b];

	/* forward search */
    for(DNAseq::const_iterator it = pattern.begin() + 1; s > 0 && it < pattern.end() && *it != 0; ++it)
    	fwdExt(p, q, s, *it);

    for(saidx_t i = p; i < p + s; ++i) { // locate fwd locs
    	saidx_t SAstart = accessSA(i);
    	locs.push_back(GLoc(SAstart, SAstart + pattern.length(), -1, GLoc::FWD));
    }
    return locs;
}

vector<GLoc> FMDIndex::locateAllRev(const DNAseq& pattern) const {
	vector<GLoc> locs;
	if(pattern.empty())
		return locs;

	nt16_t b = pattern.back();
	saidx_t p = C[b];
	saidx_t q = C[DNAalphabet::complement(b)];
	saidx_t s = C[b + 1] - C[b];

	/* backward search */
    for(DNAseq::const_reverse_iterator it = pattern.rbegin() + 1; s > 0 && it < pattern.rend(); ++it)
    	backExt(p, q, s, *it);

    for(saidx_t i = q; i < q + s; ++i) { // locate rev locs
    	saidx_t SAstart = accessSA(i);
    	locs.push_back(GLoc(SAstart, SAstart + pattern.length(), -1, GLoc::REV));
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

	const saidx_t N1 = lhs.length();
	const saidx_t N2 = rhs.length();
	const saidx_t N = N1 + N2;
	FMDIndex::BCarray_t BMerged;
	FMDIndex::BCarray_t CMerged;
	/* build merged B and C */
	for(saidx_t i = 0; i <= DNAalphabet::NT16_MAX; ++i) {
		BMerged[i] = lhs.B[i] + rhs.B[i];
		CMerged[i] = lhs.C[i] + rhs.C[i];
	}
	/* build interleaving bitstr */
	BitStr32 bstrM(N);
#pragma omp parallel for
	for(saidx_t i = 0; i < lhs.C[0 + 1]; ++i) { // i-th pass LF-mapping on lhs
		saidx_t j = i;
		saidx_t RA = rhs.C[0 + 1];
		sauchar_t b;
		do {
#pragma omp critical(WRITE_bstrM)
			bstrM.set(j + RA);
			b = lhs.bwt.access(j);
			/* LF-mapping */
			RA = rhs.LF(b, RA - 1);
			j = lhs.LF(b, j) - 1;
		}
		while(b != 0);
	}

	/* build merbed BWT */
	DNAseq bwtM;
	bwtM.reserve(N);
	for(saidx_t i = 0, j = 0, k = 0; k < N; ++k)
		bwtM.push_back(bstrM.test(k) ? lhs.bwt.access(i++) : rhs.bwt.access(j++));
	return FMDIndex(BMerged, CMerged, bwtM, lhs.keepSA || rhs.keepSA /* keep SA if any of the two operands keepSA */);
}

void FMDIndex::backExt(saidx_t& p, saidx_t& q, saidx_t& s, sauchar_t b) const {
	if(!DNAalphabet::isBasic(b))
		return;
	saidx_t pN; // fwd strand backExt
	BCarray_t qB, sB;
	qB.fill(0);
	sB.fill(0);
	/* calculate new p and s */
	saidx_t O = bwt.rank(b, p - 1);
	pN = C[b] + O;
	sB[b] = bwt.rank(b, p + s - 1) - O;

	/* update q if s changes */
	if(sB[b] != s) {
		sB[0] = bwt.rank(0, p + s - 1) - bwt.rank(0, p - 1);
		for(nt16_t i = b + 1; i <= DNAalphabet::NT16_MAX; ++i) { // search from b + 1
			if(B[i] > 0) // if is a basic symbol
				sB[i] = bwt.rank(i, p + s - 1) - bwt.rank(i, p - 1);
		}
		/* new range of [q', q' + s' - 1] is a subrange of original [q, q + s] */
		/* devide q + q + s */
		qB[0] = q;
		qB[DNAalphabet::NT16_MAX] = qB[0] + sB[0];
		for(nt16_t i = DNAalphabet::NT16_MAX; i > b; --i) // only need to search till b + 1
			qB[i - 1] = qB[i] + sB[i];
		q = qB[b];
	}

	/* update p and s */
	p = pN;
	s = sB[b];
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */
