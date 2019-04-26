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

const double SeedChain::DEFAULT_MAX_LOD10 = 4;

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
		vector<size_t> chainIdx; // chain index for this DFS search
		dfsSeeds(inputSeeds, i, outputChains, chainIdx, seedIdx, maxIndel);
	}
	assert(!outputChains.empty());
	return outputChains;
}

void SeedChain::dfsSeeds(const SeedList& inputSeeds, size_t i, ChainList& outputChains,
		vector<size_t>& chainIdx, vector<bool>& seedIdx, int64_t maxIndel) {
	if(seedIdx[i]) // already checked
		return;
	chainIdx.push_back(i);
	seedIdx[i] = true;
	bool isLeaf = true;
	for(size_t j = i + 1; j < inputSeeds.size(); ++j) {
		if(!seedIdx[j] && SeedPair::isCompatitable(inputSeeds[i], inputSeeds[j], maxIndel)) { // a compatitable child
			isLeaf = false;
			dfsSeeds(inputSeeds, j, outputChains, chainIdx, seedIdx, maxIndel);
		}
	}
	if(isLeaf) { // no more child seed found, leaf reached
		SeedChain chain;
		chain.reserve(chainIdx.size()); // chain size is the current call stack size
		for(size_t k : chainIdx)
			chain.push_back(inputSeeds[k]); // add stack values as the SeedChain path
		outputChains.push_back(chain); // add this chain
	}
	chainIdx.pop_back();
}

ChainList& SeedChain::filterChains(ChainList& chains, double maxLod10) {
	const size_t N = chains.size();
	if(N <= 1)
		return chains;
	/* sort chains by loglik, so bad chains near the end */
	std::sort(chains.begin(), chains.end(),
			[](const SeedChain& lhs, const SeedChain& rhs) { return lhs.loglik() < rhs.loglik(); });

	for(size_t i = N; i > 0; --i) { // search backward
		for(size_t j = i - 1; j > 0; --j) {
			if(contained(chains[i - 1], chains[j - 1]) && chains[i - 1].log10lik() - chains[j - 1].log10lik() > maxLod10) { // a redundant chain
				chains.erase(chains.begin() + i - 1);
				break;
			}
		}
	}
	return chains;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
