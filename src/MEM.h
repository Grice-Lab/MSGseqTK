/*
 * MEM.h
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#ifndef SRC_MEM_H_
#define SRC_MEM_H_

#include <vector>
#include <cmath>
#include <algorithm>
#include "Loc.h"
#include "DNAseq.h"
#include "QualStr.h"
#include "MSGseqCleanConst.h"

namespace EGriceLab {
namespace MSGseqClean {

/*
 * Maximal Exact Match betwen a DNAseq and a FM-index database
 */
using std::vector;

struct MEM {
	/** default constructor */
	MEM() = default;

	/** construct an MEM with all info */
	MEM(uint64_t from, uint64_t to, const DNAseq* seq, const QualStr* qual,
			int64_t N, const int64_t* B,
			const vector<Loc>& locs)
	: from(from), to(to), seq(seq), qual(qual), N(N), locs(locs)
	{
		std::copy(B, B + UINT8_MAX + 1, this->B);
	}

	/** construct an MEM with all info but not locs */
	MEM(uint64_t from, uint64_t to, const DNAseq* seq, const QualStr* qual,
			int64_t N, const int64_t* B)
	: from(from), to(to), seq(seq), qual(qual), N(N)
	{
		std::copy(B, B + UINT8_MAX + 1, this->B);
	}

	/** delegating construct an MEM with all info but not qual */
	MEM(uint64_t from, uint64_t to, const DNAseq* seq,
			int64_t N, const int64_t* B,
			const vector<Loc>& locs)
	: MEM(from, to, seq, nullptr, N, B, locs)
	{  }

	/** delegating construct an MEM with all info but not qual and locs */
	MEM(uint64_t from, uint64_t to, const DNAseq* seq,
			int64_t N, const int64_t* B)
	: MEM(from, to, seq, nullptr, N, B)
	{  }

	/* member methods */
	/** get length of this MEM */
	uint64_t length() const {
		return to - from;
	}

	/** test whether this MEM is empty */
	bool empty() const {
		return length() == 0;
	}

	/** reset this MEM to initial state */
	void reset() {
		from = 0;
		to = 0;
		seq = nullptr;
		qual = nullptr;
		locs.clear();
	}

	/**
	 * get loglikelihood of this MEM
	 * @param  baseFreq  base frequency used to determine the loglik, default equal base frequency
	 * @return  loglik using the observed sequence and quality
	 */
	double loglik() const;

	/** get the likelihood of this MEM */
	double liklihood() const {
		return ::exp(loglik());
	}

	/** get the E-value of observing this MEM on a known size database */
	double evalue() const {
		return N * liklihood();
	}

	/* member fields */
	int64_t from = 0; /* 0-based relative start on seq */
	int64_t to = 0;   /* 1-based relative end on seq */
	const DNAseq* seq = nullptr;    // using pointers for fewer overhead,
	const QualStr* qual = nullptr;  // assumes their existence

	int64_t N = 0; /* database size */
	int64_t B[UINT8_MAX + 1] = { 0 }; /* database base count */

	vector<Loc> locs; /* all Loc this MEM matches to w/ reversed coordinates */

};

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_MEM_H_ */
