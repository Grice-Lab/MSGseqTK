/*
 * SeedPair.h
 *
 *  Created on: Mar 26, 2019
 *      Author: zhengqi
 */

#ifndef SRC_SEEDPAIR_H_
#define SRC_SEEDPAIR_H_

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "GLoc.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::pair;
using std::istream;
using std::ostream;

class SeedPair;
typedef vector<SeedPair> SeedList;
typedef pair<SeedList, SeedList> SeedListPE;

/**
 * a SeedPair is a pair of GLocs can be used as an alignment seed between a query read and a target chrom
 */
class SeedPair {
public:
	/* constructors */
	/** default constructor */
	SeedPair() = default;

	/** construct from given values */
	SeedPair(int64_t from, int64_t start, int64_t len, int64_t tid, GLoc::STRAND strand, double logP)
	: from(from), to(from + len), start(start), end(start + len), tid(tid), strand(strand), logP(logP)
	{  }

	/* member methods */
	/* getters and setters */
	int64_t getEnd() const {
		return end;
	}

	int64_t getFrom() const {
		return from;
	}

	double loglik() const {
		return logP;
	}

	double pvalue() const {
		return std::exp(loglik());
	}

	double evalue(int64_t L = 1) const {
		return pvalue() * L;
	}

	int64_t getStart() const {
		return start;
	}

	GLoc::STRAND getStrand() const {
		return strand;
	}

	int64_t getTid() const {
		return tid;
	}

	int64_t getTo() const {
		return to;
	}

	/** get length */
	int64_t length() const {
		return to - from;
	}

	/** test whether is empty */
	bool empty() const {
		return length() == 0;
	}

	/**
	 * save this seed to binary output
	 * @override base class method
	 */
	ostream& save(ostream& out) const;

	/**
	 * load a seed from a binary input
	 * @override base class method
	 */
	istream& load(istream& in);

	/**
	 * write this seed to text output
	 * @override base class method
	 */
	ostream& write(ostream& out) const;

	/**
	 * read this seed from text output
	 * @override base class method
	 */
	istream& read(istream& in);

private:
	/* member fields */
	int64_t from = 0; // query begin
	int64_t to = 0; // query end
	int64_t start = 0; // target begin
	int64_t end = 0; // target end
	int64_t tid = -1; // target id
	GLoc::STRAND strand = GLoc::UNK; // query strand, note target strand is always FWD
	double logP = 0; /* logliklihood */

public:
	/* static methods */
	/** get total indels between two seeds
	 * @return  positive values indicates insertions, and negative values gives deletions
	 * or INT64_MAX if both is not compatitable
	 */
	static int64_t nIndel(const SeedPair& lhs, const SeedPair& rhs) {
		if(!isCompatitable(lhs, rhs))
			return INT64_MAX;
		return (rhs.from - lhs.to) - (rhs.start - lhs.end);
	}

	/** test whether two AlnSeeds are compatitable */
	static bool isCompatitable(const SeedPair& lhs, const SeedPair& rhs) {
		return lhs.to <= rhs.from && lhs.end <= rhs.start && lhs.tid == rhs.tid && (lhs.strand & rhs.strand) != 0;
	}

	/** test whether two AlnSeeds are compatitable and ordered and with not too much indels */
	static bool isCompatitable(const SeedPair& lhs, const SeedPair& rhs, int64_t maxIndel) {
		return std::abs(nIndel(lhs, rhs)) <= maxIndel;
	}

	/** get the best (min) loglik of a SeedList */
	static double bestLoglik(const SeedList& seeds);

	/** get the best (min) pvalue of a SeedList */
	static double bestPvalue(const SeedList& seeds) {
		return std::exp(bestLoglik(seeds));
	}

	/* non-member operators */
	/** formated input */
	friend istream& operator>>(istream& in, SeedPair& seed) {
		return seed.read(in);
	}

	/** formated output */
	friend ostream& operator<<(ostream& out, const SeedPair& seed) {
		return seed.write(out);
	}

	/** relational operator */
	friend bool operator<(const SeedPair& lhs, const SeedPair& rhs);
	friend bool operator==(const SeedPair& lhs, const SeedPair& rhs);
};

inline bool operator<(const SeedPair& lhs, const SeedPair& rhs) {
	return lhs.from != rhs.from ? lhs.from < rhs.from :
			lhs.to != rhs.to ? lhs.to < rhs.to :
					lhs.tid != rhs.tid ? lhs.tid < rhs.tid :
							lhs.strand != rhs.strand ? lhs.strand < rhs.strand :
									lhs.start != rhs.start ? lhs.start < rhs.start :
											lhs.end < rhs.end;
}

inline bool operator==(const SeedPair& lhs, const SeedPair& rhs) {
	return lhs.from == rhs.from && lhs.to == rhs.to &&
			lhs.tid == rhs.tid && lhs.strand == rhs.strand &&
			lhs.start == rhs.start && lhs.end == rhs.end;
}

inline bool operator!=(const SeedPair& lhs, const SeedPair& rhs) {
	return !(lhs == rhs);
}

inline bool operator<=(const SeedPair& lhs, const SeedPair& rhs) {
	return lhs < rhs || lhs == rhs;
}

inline bool operator>(const SeedPair& lhs, const SeedPair& rhs) {
	return rhs < lhs;
}

inline bool operator>=(const SeedPair& lhs, const SeedPair& rhs) {
	return !(lhs < rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

/** template specialization for customized hash function in std namespace */
namespace std {
template<>
class hash<EGriceLab::MSGseqTK::SeedPair> {
public:
  size_t operator() (const EGriceLab::MSGseqTK::SeedPair& seed) const {
	size_t res = 0;
	res ^= seed.getFrom() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= seed.getTo() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= seed.getTid() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= seed.getStrand() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= seed.getStart() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= seed.getEnd() + 0x9e3779b9 + (res << 6) + (res >> 2);
	return res;
  }
};

} /* namespace std */

#endif /* SRC_SEEDPAIR_H_ */
