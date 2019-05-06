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
 * one SeedChain roughly represent an Alignment, although not be unique
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

	double log10lik() const {
		return loglik() / std::log(10);
	}

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
			vector<size_t>& chainIdx, vector<bool>& seedIdx, int64_t maxMismatch, int64_t maxIndel);

public:
	/* static methods */
	/**
	 * construct a SeedChain using an raw ordered SeedList, given maxIndels allowed
	 * input SeedList will be ordered and grouped into chains using DFS algorithm
	 */
	static ChainList getChains(const SeedList& inputSeeds, int64_t maxMismatch, int64_t maxIndel);

	/** test whether one chain containing another */
	static bool containing(const SeedChain& lhs, const SeedChain& rhs) {
		return  lhs.getFrom() <= rhs.getFrom() && lhs.getTo() >= rhs.getTo() &&
				lhs.getStart() <= rhs.getStart() && lhs.getEnd() >= rhs.getEnd();
	}

	/** test whether one chain is contained by another */
	static bool contained(const SeedChain& lhs, const SeedChain& rhs) {
		return containing(rhs, lhs);
	}

	/** test whether two chain is overlapping */
	static bool overlap(const SeedChain& lhs, const SeedChain& rhs) {
		return lhs.getTid() == rhs.getTid() && (lhs.getStrand() & rhs.getStrand()) != 0 &&
				lhs.getStart() < rhs.getEnd() && lhs.getEnd() > rhs.getStart();
	}

	/** filter a ChainList by removing chains with significant smaller log10lik()
	 * and smaller chains that are contained in a larger chain
	 * filtered chains will be ordered by their loglik()
	 */
	static ChainList& filter(ChainList& chains);
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_SEEDCHAIN_H_ */
