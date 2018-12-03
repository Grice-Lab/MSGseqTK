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
#include <map>
#include <deque>
#include <algorithm>
#include "DNAseq.h"
#include "Genome.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::map;
using std::istream;
using std::ostream;
using std::vector;
using std::deque;

/**
 * class representing a MetaGenome basic information
 */
class MetaGenome {
public:
	typedef map<size_t, uint64_t> GENOME_SHIFTMAP; /* genome shift in this MetaGenome */
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
	const deque<Genome>& getGenomes() const {
		return genomes;
	}

	/** get a vector copy of Genomes in this MetaGenome */
	vector<Genome> getGenomeList() const {
		return vector<Genome>(genomes.begin(), genomes.end());
	}

	/**
	 * get the genome index at given location, or -1 if not found
	 * @param loc  0-based location on the concatenated MetaGenome
	 * @return  the index of the genome that covers this loc,
	 * or -1 if not found
	 */
	size_t getGenomeIndex(uint64_t loc) const;

	/**
	 * get the chromosome index of given location, or -1 if not found
	 * @param loc  0-based on the concatenated location
	 * @return  the index of chromosome that covers this loc,
	 * or -1 if not found
	 */
	size_t getChromIndex(uint64_t loc) const;

	/**
	 * get the genome at given location, or throws out_of_range exception
	 */
	const Genome& getGenomeAtLoc(uint64_t loc) const {
		return genomes.at(getGenomeIndex(loc));
	}

	/**
	 * get the shift of the i-th genome
	 * @param i  index of the genome
	 * @return  shift or total size of genomes before this
	 */
	uint64_t getGenomeShift(size_t i) const;

	/**
	 * get the shift of the i-th genome
	 * @param i  index of the genome
	 * @return  shift or total size of genomes before this
	 */
	GENOME_SHIFTMAP getGenomeShift() const;

	/**
	 * get a unique genome/chromo location id, which is the sorted chromosome first cover this region,
	 * or -1 if not found
	 */
	size_t getLocId(uint64_t loc) const;

	/** get genome by index */
	const Genome& getGenome(size_t i) const {
		return genomes[i];
	}

	/** count the number of genomes with a given ID */
	size_t countGenome(const string& genomeId) const;

	/** check whether this genome with given ID exists */
	bool hasGenome(const string& genomeId) const;

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

	/** relational operators */
	friend bool operator==(const MetaGenome& lhs, const MetaGenome& rhs);

	/* member fields */
private:
	deque<Genome> genomes;
};

inline bool operator==(const MetaGenome& lhs, const MetaGenome& rhs) {
	return lhs.genomes == rhs.genomes;
}

inline bool operator!=(const MetaGenome& lhs, const MetaGenome& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
