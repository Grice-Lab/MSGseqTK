/*
 * ScoreScheme.h
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#ifndef SRC_SCORESCHEME_H_
#define SRC_SCORESCHEME_H_

#include <cstdint>
#include <array>
#include "DNAalphabet.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::array;

/**
 * An alignment scoring scheme for DNAalphabet
 */
struct ScoreScheme {
	/* constructors */
	/** construct from all given values */
	ScoreScheme(uint32_t matchScore, uint32_t mismatchPenalty, uint32_t gapOPenalty, uint32_t gapEPenalty, uint32_t clipPenalty);

	/** default constructor */
	ScoreScheme()
	: ScoreScheme(DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_PENALTY, DEFAULT_GAP_OPEN_PENALTY, DEFAULT_GAP_EXT_PENALTY, DEFAULT_CLIP_PENALTY)
	{  }

	/* member methods */
	/** get match score between two DNA::alpahbet bases */
	uint32_t getScore(int8_t b1, int8_t b2) const {
		return SCORE[b1][b2];
	}

	/** get affine gap penalty of a given length */
	uint32_t gapPenalty(uint32_t length) const {
		if(length == 0)
			return 0;
		else
			return gapOPenalty + (gapEPenalty - 1) * length;
	}

	/** get affine gap penalty of opening a gap */
	uint32_t openGapPenalty() const {
		return gapOPenalty + gapEPenalty;
	}

	/** get affine gap penalty of extending a gap */
	uint32_t extGapPenalty() const {
		return gapEPenalty;
	}

	/* member fields */
	array<array<uint32_t, DNAalphabet::SIZE>, DNAalphabet::SIZE> SCORE; /* match/mismatch score matrix */
	uint32_t gapOPenalty;
	uint32_t gapEPenalty;
	uint32_t clipPenalty;

	/* static fields */
	static const uint32_t DEFAULT_MATCH_SCORE = 1;
	static const uint32_t DEFAULT_MISMATCH_PENALTY = 4;
	static const uint32_t DEFAULT_GAP_OPEN_PENALTY = 5;
	static const uint32_t DEFAULT_GAP_EXT_PENALTY = 3;
	static const uint32_t DEFAULT_CLIP_PENALTY = 5; /* 5' and 3' soft-clip (S) penalty for reads */
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_SCORESCHEME_H_ */
