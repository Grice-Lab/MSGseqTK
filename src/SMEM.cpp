/*
 * SMEM.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */

#include "SMEM.h"

namespace EGriceLab {
namespace MSGseqTK {

const double SMEM::DEFAULT_MAX_EVALUE = 0.01;

SMEM& SMEM::evaluate() {
	if(empty())
		return *this;
	logP = 0;
	int64_t N = 0;
	for(int64_t i = from; i < to; ++i) {
		DNAseq::value_type b = seq->getBase(i);
		if(fmdidx->getBaseCount(b) > 0) { // ignore non-existing bases
			logP += std::log(fmdidx->getBaseCount(b));
			N++;
		}
	}
	logP -= N * log(fmdidx->length()); /* subtract denominator */
	return *this;
}

int64_t SMEM::dbDist(const SMEM& lhs, const SMEM& rhs) {
	assert(lhs.mtg == rhs.mtg);
	int64_t minD = INT64_MAX;
	for(const Loc& loc1 : lhs.locs)
		for(const Loc& loc2 : rhs.locs)
			if(isCompatitable(lhs.mtg, loc1, loc2))
				minD = std::min(minD, Loc::dist(loc1, loc2));
	return minD;
}

ostream& SMEM::write(ostream& out) const {
	/* write basic info */
	out << from << '-' << to << ':' << loglik() << ':'; /* SAstart and SAend are ignored */
	/* write Loc info */
	for(vector<GLoc>::const_iterator loc = locs.begin(); loc != locs.end(); ++loc) {
		if(loc != locs.begin())
			out << ',';
		out << *loc;
	}
	return out;
}

SMEM_LIST SMEM::findSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		int64_t& from, int64_t& to) {
	const size_t L = seq->length();
	assert(from < L);
	SMEM_LIST curr, prev, match;

	int64_t i = from;
	nt16_t b = seq->getBase(i);
	if(!DNAalphabet::isBasic(b)) // first base is non-basic, no-matches
		return curr;

	/* forward extension */
	int64_t p = fmdidx->getCumCount(b);
	int64_t q = fmdidx->getCumCount(DNAalphabet::complement(b));
	int64_t s = fmdidx->getCumCount(b + 1) - fmdidx->getCumCount(b);
	int64_t p0 = p;
	int64_t q0 = q;
	int64_t s0 = s;
	for(i = from + 1; i <= L; ++i) {
		if(i == L) {
			curr.push_back(SMEM(seq, mtg, fmdidx, from, i, p0, q0, s0));
			break;
		}
		else {
			fmdidx->fwdExt(p, q, s, seq->getBase(i));
			if(s != s0) // a different [p, q, s] Bi-directional interval found
				curr.push_back(SMEM(seq, mtg, fmdidx, from, i, p0, q0, s0));
			if(s <= 0)
				break;
			// updates
			p0 = p;
			q0 = q;
			s0 = s;
		}
	}
	// set to
	to = i;

	/* backward extension */
	if(from == 0) {
		match = curr; // back-ext not possible
	}
	else {
		std::reverse(curr.begin(), curr.end()); // put larger SMEM in the front
		std::swap(curr, prev);
		size_t i0 = L;
		for(i = from - 1; i >= -1; --i) {
			curr.clear();
			int64_t s1 = -1;
			for(const SMEM& smem : prev) {
				p = smem.fwdStart;
				q = smem.revStart;
				s = smem.size;
				p0 = p;
				q0 = q;
				s0 = s;
				fmdidx->backExt(p, q, s, (i >= 0 ? seq->getBase(i) : 0));
				if((s <= 0 || i == -1) && curr.empty() && i < i0) {
					match.push_back(SMEM(seq, mtg, fmdidx, i + 1, smem.to, p0, q0, s0));
					i0 = i;
				}
				if(s > 0 && s1 != s) {
					curr.push_back(SMEM(seq, mtg, fmdidx, i, smem.to, p, q, s));
					s1 = s;
				}
			}
			if(curr.empty())
				break;
			std::reverse(curr.begin(), curr.end()); // put larger SMEM in the front
			std::swap(curr, prev);
		}
		from = i + 1; // update from
	}
	return match;
}

SMEM& SMEM::findLocs(size_t maxNLocs) {
	locs.reserve(size);
	const size_t N = std::min(maxNLocs, static_cast<size_t>(size));
	for(size_t i = 0; i < N; ++i) {
		{
			int64_t start = fmdidx->accessSA(fwdStart + i);
			if(mtg->getStrand(start) == GLoc::FWD) // always only search loc on fwd tStrand
				locs.push_back(GLoc(start, start + length(), mtg->getLocId(start), GLoc::FWD));
		}
		{
			int64_t start = fmdidx->accessSA(revStart + i);
			if(mtg->getStrand(start) == GLoc::FWD) // always only search loc on fwd tStrand
				locs.push_back(GLoc(start, start + length(), mtg->getLocId(start), GLoc::REV));
		}
	}
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
