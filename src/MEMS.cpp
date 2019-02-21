/*
 * MEMS.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: zhengqi
 */

#include <unordered_set>
#include <stack>
#include <cmath>
#include <utility>
#include "MEMS.h"
#include "MSGseqTKConst.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace EGriceLab {
namespace MSGseqTK {

using std::unordered_set;
using std::stack;

const double MEMS::DEFAULT_MAX_EVALUE = 0.01;

int64_t MEMS::length() const {
	int64_t L = 0;
	for(const MEM& mem : *this)
		L += mem.length();
	return L;
}

MEMS& MEMS::operator+=(const MEMS& other) {
	for(const MEM& mem : other)
		if(MEM::isCompatitable(back(), mem))
	insert(end(), mem);
	return *this;
}

MEMS MEMS::sampleMEMSfwd(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		RNG& rng, double maxEvalue) {
	MEMS mems;
	int64_t from = 0;
	while(from < seq->length()) {
		MEM mem = MEM::findMEMfwd(seq, mtg, fmdidx, from).evaluate();
		/* accept by chance */
		if(mem_dist(rng) <= maxEvalue / mem.evalue())
			mems.push_back(mem);
		from = mem.to + 1;
	}
	return mems;
}

MEMS MEMS::sampleMEMSrev(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
		RNG& rng, double maxEvalue) {
	MEMS mems;
	int64_t to = seq->length();
	while(to > 0) {
		MEM mem = MEM::findMEMrev(seq, mtg, fmdidx, to).evaluate();
		/* accept by chance */
		if(mem_dist(rng) <= maxEvalue / mem.evalue())
			mems.push_back(mem);
		to = mem.from - 1;
	}
	return mems;
}

size_t MEMS::bestMEMIndex() const {
	size_t bestIdx = -1;
	double maxLoglik = inf;
	for(MEMS::size_type i = 0; i < size(); ++i)
		if((*this)[i].loglik() < maxLoglik)
			bestIdx = i;
	return bestIdx;
}

ostream& MEMS::write(ostream& out) const {
	if(empty())
		return out;
	for(const MEM& mem : *this)
		out << mem << ';';
	out << "loglik=" << loglik();
	return out;
}

} /* namespace UCSC */
} /* namespace EGriceLab */
