/*
 * ParingScheme.cpp
 *
 *  Created on: Jun 19, 2019
 *      Author: zhengqi
 */

#include "PairingScheme.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {

const double PairingScheme::DEFAULT_MEAN_INSERT = 200;
const double PairingScheme::DEFAULT_CV_INSERT = 0.2;
const double PairingScheme::DEFAULT_OUTLIER_DIST = 10;

double PairingScheme::pr(double x) const {
	if(m == 0) // use uniform pairing prob
		return 1;
	else if(min <= x && x <= max) // in given range
		return pdf(pairingDist, x);
	else
		return 0; // not possible
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
