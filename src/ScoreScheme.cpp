/*
 * ScoreScheme.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include "ScoreScheme.h"

namespace EGriceLab {
namespace MSGseqTK {

const double ScoreScheme::DEFAULT_MATCH_SCORE = 1;
const double ScoreScheme::DEFAULT_MISMATCH_PENALTY = 4;
const double ScoreScheme::DEFAULT_GAP_OPEN_PENALTY = 5;
const double ScoreScheme::DEFAULT_GAP_EXT_PENALTY = 3;
const double ScoreScheme::DEFAULT_CLIP_PENALTY = 5; /* 5' and 3' soft-clip (S) penalty for reads */

ScoreScheme::ScoreScheme(double matchScore, double mismatchPenalty, double gapOPenalty, double gapEPenalty, double clipPenalty)
: gapOPenalty(gapOPenalty), gapEPenalty(gapEPenalty), clipPenalty(clipPenalty) {
	for(uint32_t i = 0; i < DNAalphabet::SIZE; ++i)
		for(uint32_t j = 0; j < DNAalphabet::SIZE; ++j)
			SCORE[i][j] = i == j ? matchScore : -mismatchPenalty;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

