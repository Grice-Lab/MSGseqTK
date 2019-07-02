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

const double SeedChain::MAX_LOD10 = 10;

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

ChainList SeedChain::getChains(const SeedList& inputSeeds, int64_t maxMismatch, int64_t maxIndel) {
	const size_t N = inputSeeds.size();
	vector<bool> seedIdx(N);
	ChainList outputChains;
	for(size_t i = 0; i < N; ++i) {
		vector<size_t> chainIdx; // chain index for this DFS search
		dfsSeeds(inputSeeds, i, outputChains, chainIdx, seedIdx, maxMismatch, maxIndel);
	}
	return outputChains;
}

void SeedChain::dfsSeeds(const SeedList& inputSeeds, size_t i, ChainList& outputChains,
		vector<size_t>& chainIdx, vector<bool>& seedIdx, int64_t maxMismatch, int64_t maxIndel) {
	if(seedIdx[i]) // already checked
		return;
	chainIdx.push_back(i);
	seedIdx[i] = true;
	bool isLeaf = true;
	for(size_t j = i + 1; j < inputSeeds.size(); ++j) {
		if(!seedIdx[j] && SeedPair::nMismatch(inputSeeds[i], inputSeeds[j]) <= maxMismatch &&
				std::abs(SeedPair::nIndel(inputSeeds[i], inputSeeds[j])) <= maxIndel) { // a compatitable child
			isLeaf = false;
			dfsSeeds(inputSeeds, j, outputChains, chainIdx, seedIdx, maxMismatch, maxIndel);
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

ChainList& SeedChain::uniq(ChainList& chains) {
	if(chains.size() <= 1)
		return chains;
	for(ChainList::iterator i = chains.end(); i > chains.begin(); --i) { // search backward
		for(ChainList::iterator j = i - 1; j > chains.begin(); --j) {
			if(contained(*(i - 1), *(j - 1))) { // a redundant chain
				chains.erase(i - 1);
				break;
			}
		}
	}
	assert(!chains.empty());
	return chains;
}

ChainList& SeedChain::filter(ChainList& chains, double maxLod0) {
	/* sort chains by loglik */
	std::sort(chains.begin(), chains.end(),
			[](const SeedChain& lhs, const SeedChain& rhs) { return lhs.loglik() < lhs.loglik(); });
	const double bestLog10lik = chains.front().log10lik();
	chains.erase(std::remove_if(chains.begin(), chains.end(),
			[=](const SeedChain& chain) { return chain.log10lik() > bestLog10lik + maxLod0; }),
			chains.end());
	return chains;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
