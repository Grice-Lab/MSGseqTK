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
#include <algorithm>
#include "MetaGenome.h"
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
	/** default constructor */
	MetaGenomeAnno() = default;

	/** construct an MetaGenomeAnno for a given MetaGenome, while each Genome has empty annotations */
	explicit MetaGenomeAnno(const MetaGenome& mtg);

	/** get total number of genomes */
	size_t numAnnotated() const {
		return genomeAnnos.size();
	}

	/** get total number of annotations */
	size_t numAnnotations() const;

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

	/** get a GenomeAnno iterator by id */
	vector<GenomeAnno>::const_iterator getAnno(const string& id) const {
		return std::find_if(genomeAnnos.begin(), genomeAnnos.end(), [&id](const GenomeAnno& anno) { return anno.getGenome().getId() == id; });
	}

	/** get a GenomeAnno iterator by id, non-const version*/
	vector<GenomeAnno>::iterator getAnno(const string& id) {
		return std::find_if(genomeAnnos.begin(), genomeAnnos.end(), [&id](GenomeAnno& anno) { return anno.getGenome().getId() == id; });
	}

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

	/** write this object to text output in GFF format */
	ostream& write(ostream& out) const;

	/** read aggregated GFF annotations for this MetaGenomeAnno */
	istream& read(istream& in);

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
	static string readAll(istream& in) {
		return string(std::istreambuf_iterator<char>(in), { });
	}

	/**
	 * write GFF header comments
	 */
	static ostream& writeGFFHeader(ostream& out, const string& dbName, GFF::Version ver = GFF::GFF3);

	/**
	 * read pre-built GFF header comments
	 */
	static istream& readGFFHeader(istream& in, string& dbName, GFF::Version& ver);
};

inline MetaGenomeAnno& MetaGenomeAnno::operator+=(const MetaGenomeAnno& other) {
	genomeAnnos.insert(genomeAnnos.end(), other.genomeAnnos.begin(), other.genomeAnnos.end());
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOMEANNO_H_ */
