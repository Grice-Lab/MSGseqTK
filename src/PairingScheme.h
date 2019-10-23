/*
 * PairingScheme.h
 *
 *  Created on: Jun 19, 2019
 *      Author: zhengqi
 */

#ifndef SRC_PAIRINGSCHEME_H_
#define SRC_PAIRINGSCHEME_H_

#include <cmath>
#include <cstdint>
#include <cassert>
#include <boost/math/distributions/normal.hpp>
#include "Loc.h"

namespace EGriceLab {
namespace MSGseqTK {

using boost::math::normal;

/**
 * A Pairing scheme for evaluating Alignment pairs
 * it uses a truncated normal distribution to calculate the pairing probabilities
 */
class PairingScheme {
public:
	/* constructors */
	/** default constructor */
	PairingScheme() : m(DEFAULT_MEAN_INSERT) {
		updateParam();
		updateRange();
	}

	/* member methods */
	/** getters and setters */
	double getMean() const {
		return m;
	}

	double getSD() const {
		return s;
	}

	double getMin() const {
		return min;
	}

	double getMax() const {
		return max;
	}

	void setMean(double m) {
		this->m = m;
	}

	void setSD(double s) {
		this->s = s;
	}

	void setMin(double min) {
		this->min = min;
	}

	/** set max insert */
	void setMax(double max) {
		this->max = max;
	}

	/** set insert range */
	void setRange(double min, double max) {
		this->min = min;
		this->max = max;
	}

	/** update distribution parameters */
	void updateParam() {
		if(m == 0)
			s = 1; // use standard normal to prevent zero sd error
		if(s == 0)
			s = m * DEFAULT_CV_INSERT;
		pairingDist = normal(m, s);
	}

	/** update range parameters */
	void updateRange() {
		if(m > 0 && s > 0) {
			min = m - s * DEFAULT_OUTLIER_DIST;
			max = m + s * DEFAULT_OUTLIER_DIST;
		}
	}

	/**
	 * get pvalue of observing an insert size
	 * return truncated pdf of Normal distribution, or 1 if m == 0
	 */
	double pr(double x) const;

	/** get loglik of observing an insert size */
	double loglik(double x) const {
		return std::log(pr(x));
	}

	/** get log10lik of observing an insert size */
	double log10lik(double x) const {
		return std::log10(pr(x));
	}

	/* member fields */
private:
	double m = 0;
	double s = 0;
	double min = 0;
	double max = 0;
	normal pairingDist; // underlying normal distribution

	/* static members */
public:
	static const double DEFAULT_MEAN_INSERT;
	static const double DEFAULT_CV_INSERT; // coefficient of variation
	static const double DEFAULT_OUTLIER_DIST; // outlier distance of sigma away from the mean for min and max insert
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_PAIRINGSCHEME_H_ */
