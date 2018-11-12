/*
 * MetaGenomeAnno.h
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */

#ifndef SRC_METAGENOMEANNO_H_
#define SRC_METAGENOMEANNO_H_

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include "GenomeAnno.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::istream;
using std::ostream;
using std::vector;

/**
 * class for MetaGenome annotations
 */
class MetaGenomeAnno {
public:
	/* constructors */

	/** get total number of genomes */
	size_t numAnnotated() const {
		return genomeAnnos.size();
	}

	/**
	 * add an annotation at the end
	 */
	void push_back(const GenomeAnno& anno) {
		genomeAnnos.push_back(anno);
	}

	/** append multiple annotations at the end */
	void append(const vector<GenomeAnno>& blockAnnos) {
		genomeAnnos.insert(genomeAnnos.end(), blockAnnos.begin(), blockAnnos.end());
	}

	/** write this object to text output in GFF format */
	ostream& write(ostream& out) const;

	/** merge this MetaGenomeAnno with another one,
	 * with its name unchanged
	 */
	MetaGenomeAnno& operator+=(const MetaGenomeAnno& other);

	/* non-member operators */
	/** concate two MetaGenomes and return a new copy */
	friend MetaGenomeAnno operator+(const MetaGenomeAnno& lhs, const MetaGenomeAnno& rhs) {
		MetaGenomeAnno annoMerged(lhs);
		annoMerged += rhs;
		return annoMerged;
	}

	/* member fields */
private:
	vector<GenomeAnno> genomeAnnos;

	/* static methods */
public:
	/** read in all content from a GFF text input */
	static string read(istream& in) {
		return string(std::istreambuf_iterator<char>(in), { });
	}
};

inline MetaGenomeAnno& MetaGenomeAnno::operator+=(const MetaGenomeAnno& other) {
	genomeAnnos.insert(genomeAnnos.end(), other.genomeAnnos.begin(), other.genomeAnnos.end());
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOMEANNO_H_ */
