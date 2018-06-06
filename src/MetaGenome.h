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
#include <deque>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include "DNAseq.h"
#include "Genome.h"

namespace EGriceLab {
namespace MSGseqClean {

using std::string;
using std::istream;
using std::ostream;
using std::vector;
using std::deque;
using Eigen::Vector4d;

/**
 * class representing a genome basic information
 */
class MetaGenome {
public:
	/* constructors */
	/** default constructor */
//	MetaGenome() = default;

	/** get total size of this MetaGenome */
	uint64_t getSize() const;

	/** get the aggregate base count */
	Vector4l getBaseCount() const;

	/** get total number of genomes */
	size_t numGenomes() const {
		return genomes.size();
	}

	/** get all Genomes in this MetaGenome */
	vector<Genome> getGenomes() const;

	/** get all genome names */
	vector<string> getGenomeNames() const;

	/**
	 * get the genome index at given location, or -1 if not found
	 * @param loc  0-based location on the MetaGenome
	 * @return  the index of the genome that covers this loc
	 */
	size_t getGenomeIndex(uint64_t loc) const;

	/**
	 * get the genome at given location, or throws out_of_range exception
	 */
	const Genome& getGenome(uint64_t loc) const {
		return genomes.at(getGenomeIndex(loc));
	}

	/**
	 * add a genome at the end of this MetaGenome
	 */
	void push_back(const Genome& genome) {
		genomes.push_back(genome);
	}

	/**
	 * add a genome at the beginning of this MetaGenome
	 */
	void push_front(const Genome& genome) {
		genomes.push_front(genome);
	}

	/** return the base-frequency array of this genome */
	Eigen::Vector4d getBaseFreq() const {
		Vector4l count = getBaseCount();
		return count.cast<double>() / static_cast<double> (count.sum());
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
	deque<Genome> genomes;
};

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
