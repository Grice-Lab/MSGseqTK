/*
 * Genome.h
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#ifndef SRC_GENOME_H_
#define SRC_GENOME_H_

#include <string>
#include <iostream>
#include <cstdint> // C++11
#include <Eigen/Dense>
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGseqClean {

using std::string;
using std::istream;
using std::ostream;
using Eigen::Vector4d;

typedef Eigen::Matrix<uint64_t, 4, 1> Vector4l;

/**
 * Basic information of a genome (within a MetaGenome)
 */
class Genome {
public:
	/* constructors */
	/** default constructor */
	Genome() = default;

	/** construct from a name and given DNAseq */
	Genome(const string& name, const DNAseq& seq);

	/* member methods */
	/* non-member functions */
	friend bool operator==(const Genome& lhs, const Genome& rhs);

	const Vector4l& getBaseCount() const {
		return baseCount;
	}

	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	uint64_t getSize() const {
		return size;
	}

	/** return the base-frequency of this genome */
	Vector4d getBaseFreq() const {
		return baseCount.cast<double>() / static_cast<double> (baseCount.sum());
	}

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

private:
	string name;
	uint64_t size = 0;
	Vector4l baseCount = Vector4l::Zero();
};

inline bool operator==(const Genome& lhs, const Genome& rhs) {
	return lhs.name == rhs.name && lhs.size == rhs.size && lhs.baseCount == rhs.baseCount;
}

inline bool operator!=(const Genome& lhs, const Genome& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_GENOME_H_ */
