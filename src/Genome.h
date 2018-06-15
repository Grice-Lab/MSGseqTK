/*
 * Genome.h
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#ifndef SRC_GENOME_H_
#define SRC_GENOME_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cstdint> // C++11
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGseqClean {

using std::string;
using std::map;
using std::vector;
using std::istream;
using std::ostream;

/**
 * Basic information of a genome (within a MetaGenome)
 */
class Genome {
	typedef map<string, uint64_t> chrmap_t;

public:
	/* constructors */
	/** default constructor */
	Genome() = default;

	/** construct using a single chrom and seq */
	Genome(const string& chr, const DNAseq& seq) {
		addChrom(chr, seq);
	}

	/* member methods */
	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	/** get number of chromosomes */
	size_t numChroms() const {
		return chromSize.size();
	}

	/** get all chromosomes */
	vector<string> getChroms() const;

	/** get chromosome size map */
	const chrmap_t& getChromSize() const {
		return chromSize;
	}

	/** get a single chromsome size */
	uint64_t getChromSize(const string& chr) const {
		return chromSize.at(chr);
	}

	/** get the overall size of this genome */
	uint64_t getSize() const;

	/** add a new chromosome with given seq */
	void addChrom(const string& chr, const DNAseq& seq);

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

	/* non-member functions */
	friend bool operator==(const Genome& lhs, const Genome& rhs);

private:
	string name;
	chrmap_t chromSize;
};

inline void Genome::addChrom(const string& chr, const DNAseq& seq) {
	chromSize[chr] = seq.length();
}

inline bool operator!=(const Genome& lhs, const Genome& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_GENOME_H_ */
