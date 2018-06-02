/*
 * MetaGenome.h
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#ifndef SRC_METAGENOME_H_
#define SRC_METAGENOME_H_

#include <cstdint>
#include <string>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGseqClean {

using std::string;
using std::array;
using std::istream;
using std::ostream;
using std::vector;
using Eigen::Vector4d;

typedef Eigen::Matrix<uint64_t, 4, 1> Vector4l;

/**
 * class representing a genome basic information
 */
class MetaGenome {
public:
	/* constructors */
	/** default constructor */
//	MetaGenome() = default;

	uint64_t getSize() const {
		return size;
	}

	const Vector4l& getBaseCount() const {
		return baseCount;
	}

	const vector<string>& getGenomeNames() const {
		return genomeNames;
	}

	const string& genomeName(size_t i) const {
		return genomeNames[i];
	}

	size_t numGenomes() const {
		return genomeNames.size();
	}

	MetaGenome& addGenome(const string& genomeName, const DNAseq& genomeSeq);

	/** return the base-frequency array of this genome */
	Eigen::Vector4d getBaseFreq() const {
		return baseCount.cast<double>() / static_cast<double> (baseCount.sum());
	}

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

	/** merge this MetaGenome with another one,
	 * with its name unchanged
	 */
	MetaGenome& operator+=(const MetaGenome& other);

	/* non-member operators */
	/** concate two MetaGenomes and return a new copy */
	friend MetaGenome operator+(const MetaGenome& lhs, const MetaGenome& rhs) {
		MetaGenome mgMerged(lhs);
		mgMerged += rhs;
		return mgMerged;
	}

	/* member fields */
private:
	uint64_t size = 0;
	vector<string> genomeNames; /* individual genome names */
	Vector4l baseCount = Vector4l::Zero();
};

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
