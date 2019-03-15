/*
 * ScoreScheme.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include "ScoreScheme.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace MSGseqTK {

const double ScoreScheme::DEFAULT_MATCH_SCORE = 1;
const double ScoreScheme::DEFAULT_MISMATCH_PENALTY = 3;
const double ScoreScheme::DEFAULT_GAP_OPEN_PENALTY = 5;
const double ScoreScheme::DEFAULT_GAP_EXT_PENALTY = 2;
const double ScoreScheme::DEFAULT_CLIP_PENALTY = 4; /* 5' and 3' soft-clip (S) penalty for reads */

void ScoreScheme::updateScores() {
	for(nt16_t i = 0; i < DNAalphabet::SIZE; ++i)
		for(nt16_t j = 0; j < DNAalphabet::SIZE; ++j) {
			uint32_t nOverlap = libSDS::popcount8(i & j);
			if(nOverlap > 0) // match
				SCORE[i][j] = matchScore * nOverlap / std::max(libSDS::popcount8(i), libSDS::popcount8(j));
			else
				SCORE[i][j] = -mismatchPenalty; /* intesecting bits is a match */
		}
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

