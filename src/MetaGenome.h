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
namespace MSGseqTK {

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
	uint64_t size() const;

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
	 * get the chromosome index of given location, or -1 if not found
	 * @param loc  0-based location
	 * @return  the index of chromosome that covers this loc
	 */
	size_t getChromIndex(uint64_t loc) const;

	/**
	 * get the genome at given location, or throws out_of_range exception
	 */
	const Genome& getGenomeAtLoc(uint64_t loc) const {
		return genomes.at(getGenomeIndex(loc));
	}

	/** get genome by index */
	const Genome& getGenome(size_t i) const {
		return genomes[i];
	}

	/** count the number of genomes with a given name */
	size_t countGenome(const string& genomeName) const;

	/** check whether this genome with given name exists */
	bool hasGenome(const string& genomeName) const;

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

	/** append multiple genomes at the end of this MetaGenome */
	void append(const vector<Genome>& blockGenomes) {
		genomes.insert(genomes.end(), blockGenomes.begin(), blockGenomes.end());
	}

	/** prepend multiple genomes at the beginning of this MetaGenome */
	void prepend(const vector<Genome>& blockGenomes) {
		genomes.insert(genomes.begin(), blockGenomes.begin(), blockGenomes.end());
	}

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

	/** write this object to text output in GFF format */
	ostream& writeGFF(ostream& out, UCSC::GFF::Version ver = UCSC::GFF::GFF3, const string& src = ".") const;

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

	/**
	 * test whether two MeteGenomes are equal
	 * @return  true if and only if all their Genomes are equal and in same order
	 */
	friend bool operator==(const MetaGenome& lhs, const MetaGenome& rhs);

	/* member fields */
private:
	deque<Genome> genomes;
};

inline bool operator!=(const MetaGenome& lhs, const MetaGenome& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
