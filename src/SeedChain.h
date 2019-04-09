/*
 * SeedChain.h
 *
 *  Created on: Mar 29, 2019
 *      Author: zhengqi
 */

#ifndef SRC_SEEDCHAIN_H_
#define SRC_SEEDCHAIN_H_

#include <vector>
#include <stack>
#include <utility>
#include "SeedPair.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::stack;
using std::pair;
class SeedChain;

typedef vector<SeedChain> ChainList;
typedef pair<ChainList, ChainList> ChainListPE;

/**
 * a SeedChain IS A SeedList that is colinearly ordered and compatitable
 */
class SeedChain : public SeedList {
public:
	/* constructors */
	/** default constructor */
	SeedChain() = default;

	/* member methods */
	int64_t getFrom() const {
		return front().getFrom();
	}

	int64_t getTo() const {
		return back().getTo();
	}

	int64_t getStart() const {
		return front().getStart();
	}

	int64_t getEnd() const {
		return back().getEnd();
	}

	/** get tid */
	int64_t getTid() const {
		return front().getTid();
	}

	/** get strand */
	GLoc::STRAND getStrand() const {
		return front().getStrand();
	}

	/** get overall loglik */
	double loglik() const;

	/** get overal length */
	int64_t length() const;

	/** get overal pvalue */
	double pvalue() const {
		return std::exp(loglik());
	}

private:
	/* private util methods */
	/**
	 * a Depth-First-Search algorithm of finding all leaf paths (chains) from a raw SeedList
	 * it recursively calls itself on newly unvisited seed, or return and append to the output list
	 */
	static void dfsSeeds(const SeedList& inputSeeds, size_t i, ChainList& outputChains,
			vector<size_t>& chainIdx, vector<bool>& seedIdx, int64_t maxIndel);

public:
	/* static methods */
	/**
	 * construct a SeedChain using an raw ordered SeedList, given maxIndels allowed
	 * input SeedList will be ordered and grouped into chains using DFS algorithm
	 */
	static ChainList getChains(const SeedList& inputSeeds, int64_t maxIndel);
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_SEEDCHAIN_H_ */
