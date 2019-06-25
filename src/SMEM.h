/*
 * SSMEM.h
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */

#ifndef SSMEM_H_
#define SSMEM_H_

#include <vector>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <utility>
#include <algorithm>
#include <set>
#include "SeedPair.h"
#include "PrimarySeq.h"
#include "FMDIndex.h"
#include "MetaGenome.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::vector;
using std::pair;
using std::set;

class SMEM;
typedef SMEM MEM;

/**
 * MEM: Maximum Exact Match between a PrimarySeq read and a MetaGenome/FMD-Iidex target that cannot extend to either end
 * SMEM: Super MEM that is not contained in any other MEM
 */
class SMEM {
public:
	/** default constructor */
	SMEM() = default;

	/** construct an SMEM with all info */
	SMEM(const PrimarySeq* seq,  const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t from, int64_t to, int64_t fwdStart, int64_t revStart, int64_t size)
	: seq(seq), mtg(mtg), fmdidx(fmdidx),
	  from(from), to(to), fwdStart(fwdStart), revStart(revStart), size(size)
	{
		evaluate();
	}

	/* member methods */
	/** getters and setters */
	const PrimarySeq* getSeq() const {
		return seq;
	}

	const FMDIndex* getFmdidx() const {
		return fmdidx;
	}

	const MetaGenome* getMtg() const {
		return mtg;
	}

	int64_t getFrom() const {
		return from;
	}

	int64_t getTo() const {
		return to;
	}

	int64_t getSize() const {
		return size;
	}

	int64_t getFwdStart() const {
		return getFwdStart();
	}

	int64_t getRevStart() const {
		return getRevStart();
	}

	double loglik() const {
		return logP;
	}

	double loglik(int64_t from, int64_t to) const;

	double pvalue() const {
		return std::exp(loglik());
	}

	double pvalue(int64_t from, int64_t to) const {
		return std::exp(loglik(from, to));
	}

	double evalue() const {
		return fmdidx->length() * pvalue();
	}

	double evalue(int64_t from, int64_t to) const {
		return fmdidx->length() * pvalue(from, to);
	}

	/** get length of this SMEM */
	int64_t length() const {
		return to - from;
	}

	/** test whether this SMEM is empty */
	bool empty() const {
		return length() == 0;
	}

	/** test whether this SMEM is valid */
	bool isValid() const {
		return size > 0;
	}

	/**
	 * get SeedPairs of this SMEM
	 * @return  mapped seeds that always on FWD strand of target
	 */
	SeedList getSeeds(int64_t maxNSeed = MAX_NSEED) const;

	/**
	 * forward extend this SMEM at current end
	 * forward extension can be one pass the seq end, which guarenteed a null base
	 * @return  updated SMEM
	 */
	SMEM& fwdExt() {
		fmdidx->fwdExt(fwdStart, revStart, size, seq->getBase(to));
		fwdEvaluate();
		to++;
		return *this;
	}

	/**
	 * get a copy of forward extension of this SMEM
	 * @return  new SMEM with extended values
	 */
	SMEM fwdExt() const {
		SMEM smemExt(*this);
		return smemExt.fwdExt();
	}

	/**
	 * backward extend this SMEM at current start
	 * backward extension can be one pass seq start, which will be a null
	 * @return  updated SMEM
	 */
	SMEM& backExt() {
		fmdidx->backExt(fwdStart, revStart, size, from > 0 ? seq->getBase(from - 1) : 0);
		backEvaluate();
		from--;
		return *this;
	}

	/**
	 * get a copy of backward extension of this SMEM
	 * @return  new SMEM with extended values
	 */
	SMEM backExt() const {
		SMEM smemExt(*this);
		return smemExt.backExt();
	}

	/**
	 * evaluate the log-probality of this SMEM
	 * @return  log-pvalue of observing this SMEM by chance, using base-frequency only
	 */
	SMEM& evaluate();

	/** forward evaluation during forward extension */
	SMEM& fwdEvaluate() {
		assert(to <= seq->length());
		logP += fmdidx->loglik(seq->getBase(to));
		return *this;
	}

	/** backward evaluation during backward extension */
	SMEM& backEvaluate() {
		assert(from >= 0);
		logP += from > 0 ? fmdidx->loglik(seq->getBase(from - 1)) : fmdidx->loglik(DNAalphabet::GAP_BASE);
		return *this;
	}

	/* static methods */
	/**
	 * find list of MEMS of a given seq starting at given position using forward/backward extension
	 */
	static MEM findMEM(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t& from, int64_t& to);

private:
	/* member fields */
	const PrimarySeq* seq = nullptr;
	const MetaGenome* mtg = nullptr;
	const FMDIndex* fmdidx = nullptr;
	int64_t from = 0; /* 0-based relative start on seq */
	int64_t to = 0;   /* 1-based relative end on seq */
	int64_t fwdStart = 0; /* 0-based start position on SA of fwd match */
	int64_t revStart = 0;   /* 0-based end position on SA of rev match */
	int64_t size = 0; /* size of match region on SA */
	double logP = 0; /* log-probability (loglik) of observing this SMEM by chance */

public:
	static const int64_t MAX_NSEED = 500; // maximum size (# of occurrence of a good SMEM)
	/* non-member functions */
	/** formated output */
	friend ostream& operator<<(ostream& out, const SMEM& smem);

	/** relationship operators */
	friend bool operator<(const SMEM& lhs, const SMEM& rhs);
	friend bool operator==(const SMEM& lhs, const SMEM& rhs);

	friend class SMEM_LIST;
};

class SMEM_LIST;
typedef SMEM_LIST MEM_LIST;
typedef pair<MEM_LIST, MEM_LIST> MEM_LIST_PE;
typedef pair<SMEM_LIST, SMEM_LIST> SMEM_LIST_PE;

/**
 * a SMEM_LIST is a list SMEM with additional methods
 */
class SMEM_LIST : public vector<SMEM> {
public:
	/* member methods */
	/** get total loglik */
	double loglik() const;

	/** get total pvalue */
	double pvalue() const {
		return std::exp(loglik());
	}

	/** get total evalue */
	double evalue() const {
		return front().fmdidx->length() * pvalue();
	}

	/** get best MEM (smallest loglik) */
	const SMEM& getBest() const {
		assert(!empty());
		return *std::min_element(begin(), end(),
				[](const SMEM& lhs, const SMEM& rhs) { return lhs.loglik() < rhs.loglik(); });
	}

	/** get worst MEM (largest loglik) */
	const SMEM& getWorst() const {
		assert(!empty());
		return *std::max_element(begin(), end(),
				[](const SMEM& lhs, const SMEM& rhs) { return lhs.loglik() < rhs.loglik(); });
	}

	/** get longest MEM */
	const SMEM& getLongest() const {
		assert(!empty());
		return *std::max_element(begin(), end(),
				[](const SMEM& lhs, const SMEM& rhs) { return lhs.length() < rhs.length(); });
	}

	/** get shortest MEM */
	const SMEM& getShortest() const {
		assert(!empty());
		return *std::min_element(begin(), end(),
				[](const SMEM& lhs, const SMEM& rhs) { return lhs.length() < rhs.length(); });
	}

	/**
	 * filter this SMEMS list by evalue and size
	 */
	SMEM_LIST& filter(int64_t minLen = MIN_LENGTH) {
		erase(std::remove_if(begin(), end(),
				[=](const SMEM& mem) { return !(minLen <= mem.length() ); }),
				end());
		return *this;
	}

	/**
	 * sort this SMEM_LIST
	 */
	SMEM_LIST& sort() {
		std::sort(begin(), end());
		return *this;
	}

	/** get unique SMEM of this list */
	SMEM_LIST& uniq() {
		erase(std::unique(begin(), end()), end());
		return *this;
	}

	/** get from of an SMEM_LIST */
	int64_t getFrom() const;

	/** get to of an SMEM_LIST */
	int64_t getTo() const;

	/** get number of valid in this SMEM list */
	size_t numValid() const {
		return std::count_if(begin(), end(),
				[](const SMEM& smem) { return smem.isValid(); });
	}

	/** test whether all of SMEM is valid */
	bool allValid() const {
		return std::all_of(begin(), end(),
				[](const SMEM& smem) { return smem.isValid(); });
	}

	/** test whether any of SMEM is valid */
	bool anyValid() const {
		return std::any_of(begin(), end(),
				[](const SMEM& smem) { return smem.isValid(); });
	}

	/** test whether none of SMEM is valid */
	bool noneValid() const {
		return std::none_of(begin(), end(),
				[](const SMEM& smem) { return smem.isValid(); });
	}

	/* static methods */
	/**
	 * find longest SMEMS of a given seq using step-wise forward/backward searches
	 */
	static SMEM_LIST findMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t minLen = MIN_LENGTH, double maxEvalue = MAX_EVALUE);

	/**
	 * find all SMEMS of a given seq starting at given position relative to the seq by forward/backward extensions
	 * @return  a SMEM_LIST overlappling position from that may contain duplicated copies
	 */
	static SMEM_LIST findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t& from, int64_t& to,
			int64_t minLen, double maxEvalue = MAX_EVALUE, int64_t minSize = MIN_SIZE);

public:
	/**
	 * find SMEMS of a given seq by step-wise searches
	 * @return SMEMS sorted by their loglik() smaller loglik near the beginning
	 */
	static SMEM_LIST findAllSMEMS(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t minLen = MIN_LENGTH, int64_t maxLen = MAX_LENGTH, double maxEvalue = MAX_EVALUE);

	/**
	 * get a SeedList of a given seq using step-wise forward/backward searches
	 * seeds will be filtered and sorted
	 */
	static SeedList findSeeds(const PrimarySeq* seq, const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t minLen = MIN_LENGTH, int64_t maxLen = MAX_LENGTH,
			double maxEvalue = MAX_EVALUE, int64_t maxNSeed = SMEM::MAX_NSEED, bool discardSeed = false);

	/**
	 * find MEMS_PE for paired-end reads
	 */
	static MEM_LIST_PE findMEMS_PE(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq,
			const MetaGenome* mtg, const FMDIndex* fmdidx,
			int64_t minLen = MIN_LENGTH, double maxEvalue = MAX_EVALUE) {
		return MEM_LIST_PE(findMEMS(fwdSeq, mtg, fmdidx, minLen, maxEvalue),
				findMEMS(revSeq, mtg, fmdidx, minLen, maxEvalue));
	}

	/**
	 * find SeedListPE for pair-end reads
	 */
	static SeedListPE findSeedsPE(const PrimarySeq* fwdSeq, const PrimarySeq* revSeq,
			 const MetaGenome* mtg, const FMDIndex* fmdidx,
			 int64_t minLen = MIN_LENGTH, int64_t maxLen = MAX_LENGTH, double maxEvalue = MAX_EVALUE,
			 int64_t maxNSeed = SMEM::MAX_NSEED, bool discardSeed = false) {
		return SeedListPE(findSeeds(fwdSeq, mtg, fmdidx, minLen, maxLen, maxEvalue, maxNSeed, discardSeed),
				findSeeds(revSeq, mtg, fmdidx, minLen, maxLen, maxEvalue, maxNSeed, discardSeed));
	}

	/** get loglik for SMEMS_PE */
	static double loglik(const SMEM_LIST_PE& smemsPE) {
		return smemsPE.first.loglik() + smemsPE.second.loglik();
	}

	/** append a new SMEM_LIST */
	SMEM_LIST& operator+=(const SMEM_LIST& other) {
		insert(end(), other.begin(), other.end());
		return *this;
	}

	/** non-member operators */
	friend SMEM_LIST operator+(const SMEM_LIST& lhs, const SMEM_LIST& rhs);

	/* static fields */
	static const int64_t MIN_LENGTH = 17; // minimum length for a significant SMEM
	static const int64_t MAX_LENGTH = 0; // maximum length for a significant SMEM
	static const int64_t MIN_SIZE = 1; // minimum size (# of occurrence of a significant SMEM)
	static const double MAX_EVALUE;
};

inline ostream& operator<<(ostream& out, const SMEM& smem) {
	out << smem.from << "-" << smem.to << ":" << smem.size << ":" << smem.fwdStart << ":" << smem.revStart << ":" << smem.logP;
	return out;
}

inline bool operator<(const SMEM& lhs, const SMEM& rhs) {
	assert(lhs.seq == rhs.seq && lhs.mtg == rhs.mtg && lhs.fmdidx == rhs.fmdidx);
	/* from and to determine all other values including size, fwdStart, revStart */
	return lhs.from != rhs.from ? lhs.from < rhs.from : lhs.to < rhs.to;
}

inline bool operator==(const SMEM& lhs, const SMEM& rhs) {
	assert(lhs.seq == rhs.seq && lhs.mtg == rhs.mtg && lhs.fmdidx == rhs.fmdidx);
	return lhs.from == rhs.from && lhs.to == rhs.to;
}

inline bool operator!=(const SMEM& lhs, const SMEM& rhs) {
	return !(lhs == rhs);
}

inline bool operator<=(const SMEM& lhs, const SMEM& rhs) {
	return lhs < rhs || lhs == rhs;
}

inline bool operator>(const SMEM& lhs, const SMEM& rhs) {
	return rhs < lhs;
}

inline bool operator>=(const SMEM& lhs, const SMEM& rhs) {
	return !(lhs < rhs);
}

inline SMEM_LIST operator+(const SMEM_LIST& lhs, const SMEM_LIST& rhs) {
	SMEM_LIST smems(lhs);
	smems += rhs;
	return smems;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

namespace std {

/**
 * template specialization for customized hash function in std namespace
 * the hash code for SMEM is only based by its pointers and from and to, where all other fields are determined
 */
template<>
class hash<EGriceLab::MSGseqTK::SMEM> {
public:
  size_t operator() (const EGriceLab::MSGseqTK::SMEM& smem) const {
	size_t res = 0;
	res ^= reinterpret_cast<uintptr_t>(smem.getSeq()) + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= reinterpret_cast<uintptr_t>(smem.getMtg()) + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= reinterpret_cast<uintptr_t>(smem.getFmdidx()) + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= smem.getFrom() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= smem.getTo() + 0x9e3779b9 + (res << 6) + (res >> 2);
//	res ^= smem.getSize() + 0x9e3779b9 + (res << 6) + (res >> 2);
//	res ^= smem.getFwdStart() + 0x9e3779b9 + (res << 6) + (res >> 2);
//	res ^= smem.getRevStart() + 0x9e3779b9 + (res << 6) + (res >> 2);
	return res;
  }
};

} /* namespace std */

#endif /* SSMEM_H_ */
