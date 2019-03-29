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
class ScoreScheme {
public:
	/* constructors */
	/** default constructor */
	ScoreScheme() {
		updateScores();
	}

	/** construct from all given values */
	ScoreScheme(double matchScore, double mismatchPenalty, double gapOPenalty, double gapEPenalty, double clipPenalty)
	: matchScore(matchScore), mismatchPenalty(mismatchPenalty),
	  gapOPenalty(gapOPenalty), gapEPenalty(gapEPenalty), clipPenalty(clipPenalty)
	{
		updateScores();
	}

	/* member methods */


	/* update match/mismatch scores
	 * matched score is weighted by # of overlapping bits to relative down-weight IUPAC ambiguous bases */
	void updateScores();

	/** get match score between two DNA::alpahbet bases */
	double getScore(int8_t b1, int8_t b2) const {
		return SCORE[b1][b2];
	}

	/** get affine gap penalty of a given length */
	double gapPenalty(uint32_t length) const {
		if(length == 0)
			return 0;
		else
			return gapOPenalty + (gapEPenalty - 1) * length;
	}

	/** get affine gap penalty of opening a gap */
	double openGapPenalty() const {
		return gapOPenalty + gapEPenalty;
	}

	/** get affine gap penalty of extending a gap */
	double extGapPenalty() const {
		return gapEPenalty;
	}

	/** set match score */
	void setMatchScore(double matchScore) {
		this->matchScore = matchScore;
		updateScores();
	}

	/** set mis-match penalty */
	void setMismatchPenalty(double mismatchPenalty) {
		this->mismatchPenalty = mismatchPenalty;
		updateScores();
	}

	/* member fields */
private:
	array<array<double, DNAalphabet::SIZE>, DNAalphabet::SIZE> SCORE; /* match/mismatch score matrix */
	double matchScore = DEFAULT_MATCH_SCORE;
	double mismatchPenalty = DEFAULT_MISMATCH_PENALTY;
	double gapOPenalty = DEFAULT_GAP_OPEN_PENALTY;
	double gapEPenalty = DEFAULT_GAP_EXT_PENALTY;
	double clipPenalty = DEFAULT_CLIP_PENALTY;

public:
	/* static fields */
	static const double DEFAULT_MATCH_SCORE;
	static const double DEFAULT_MISMATCH_PENALTY;
	static const double DEFAULT_GAP_OPEN_PENALTY;
	static const double DEFAULT_GAP_EXT_PENALTY;
	static const double DEFAULT_CLIP_PENALTY; /* 5' and 3' soft-clip (S) penalty for reads */
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_SCORESCHEME_H_ */
