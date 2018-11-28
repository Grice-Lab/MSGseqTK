/*
 * ScoreScheme.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include "ScoreScheme.h"

namespace EGriceLab {
namespace MSGseqTK {

ScoreScheme::ScoreScheme(uint32_t matchScore, uint32_t mismatchPenalty, uint32_t gapOPenalty, uint32_t gapEPenalty, uint32_t clipPenalty)
: gapOPenalty(gapOPenalty), gapEPenalty(gapEPenalty), clipPenalty(clipPenalty) {
	for(uint32_t i = 0; i < DNAalphabet::SIZE; ++i)
		for(uint32_t j = 0; j < DNAalphabet::SIZE; ++j)
			SCORE[i][j] = i == j ? matchScore : mismatchPenalty;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

