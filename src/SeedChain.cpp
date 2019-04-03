/*
 * SeedChain.cpp
 *
 *  Created on: Mar 29, 2019
 *      Author: zhengqi
 */

#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include "SeedChain.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::stack;

double SeedChain::loglik() const {
	double ll = 0;
	for(const SeedPair& seed : *this)
		ll += seed.loglik();
	return ll;
}

int64_t SeedChain::length() const {
	int64_t len = 0;
	for(const SeedPair& seed : *this)
		len += seed.length();
	return len;
}

ChainList SeedChain::getChains(const SeedList& inputSeeds, int64_t maxIndel) {
	assert(!inputSeeds.empty());
	const size_t N = inputSeeds.size();
	vector<bool> seedIdx(N);
	ChainList outputChains;
	for(size_t i = 0; i < N; ++i) {
		if(!seedIdx[i]) { // new root
			vector<size_t> chainIdx; // chain index for this DFS search
			chainIdx.push_back(i); // add root
			dfsSeeds(inputSeeds, outputChains, chainIdx, seedIdx, maxIndel);
		}
	}
	assert(!outputChains.empty());
	return outputChains;
}

void SeedChain::dfsSeeds(const SeedList& inputSeeds, ChainList& outputChains,
		vector<size_t>& chainIdx, vector<bool>& seedIdx, int64_t maxIndel) {
	assert(!chainIdx.empty());
	const size_t N = inputSeeds.size();
	size_t i = chainIdx.back();
	for(size_t j = i + 1; j < N; ++j) {
		if(!seedIdx[j] && SeedPair::isCompatitable(inputSeeds[i], inputSeeds[j], maxIndel)) { // a compatitable child
			seedIdx[j] = true;
			chainIdx.push_back(j);
			dfsSeeds(inputSeeds, outputChains, chainIdx, seedIdx, maxIndel);
		}
	}
	/* returning from recursive calls, add current stack elements in order */
	const size_t n = outputChains.size();
	for(size_t i : chainIdx)
		outputChains[n].push_back(inputSeeds[i]);
	chainIdx.pop_back();
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
