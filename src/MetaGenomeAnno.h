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
#include "GFF.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::istream;
using std::ostream;
using std::vector;
using UCSC::GFF;

/**
 * class for MetaGenome annotations
 */
class MetaGenomeAnno {
public:
	typedef map<string, vector<GFF>> GENOME_ANNOMAP; // per-genome annos
	typedef map<string, vector<GFF>> CHROM_ANNOMAP;  // per-chrom annos
	/* constructors */
	/** default constructor */
	MetaGenomeAnno() = default;

	/** construct MetaGenomeAnno from MetaGenome, and add genome-level annos */
	explicit MetaGenomeAnno(const MetaGenome& mtg);

	/** get total number of annotated genomes */
	size_t numAnnotatedGenomes() const {
		return genomeAnnos.size();
	}

	/** get total number of annotated chromosomes */
	size_t numAnnotatedChroms() const {
		return chromAnnos.size();
	}

	/** get total number of annotations */
	size_t numAnnotations() const;

	/** add a Genome to this annotation, create associated per-genome annotations */
	void addGenome(const Genome& genome);

	/** add a new genome-level GFF annotation */
	void addGenomeAnno(const string& id, const GFF& gff) {
		genomeAnnos[id].push_back(gff);
	}

	/**
	 * add a new chrom-level GFF annotations
	 * @param chrId  a unique id for a chrom
	 */
	void addChromAnno(const string& chrId, const GFF& gff) {
		chromAnnos[chrId].push_back(gff);
	}

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

	/** write annotations to GFF text output for given genome */
	ostream& write(ostream& out, const Genome& genome) const;

	/** read annotations and associated genome from GFF text input for given genomeId */
	istream& read(istream& in, const Genome& gnome, GFF::Version ver = ANNO_VER);

	/** write all annotations to GFF text output */
	ostream& write(ostream& out, const MetaGenome& mtg) const;

	/** read all GFF annotations from GFF text input */
	istream& read(istream& in, const MetaGenome& mtg);

	/**
	 * merge this MetaGenomeAnno with another one,
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

	/** relational operators */
	friend bool operator==(const MetaGenomeAnno& lhs, const MetaGenomeAnno& rhs);

	/* member fields */
private:
	GENOME_ANNOMAP genomeAnnos; // genome id->vector<GFF>
	CHROM_ANNOMAP chromAnnos;   // chrom id (genomeId.chrName)->vector<GFF>

public:
	/* static member fields */
	static const GFF::Version ANNO_VER = GFF::GFF3;
	static const string GENOME_START_TAG;
	static const string GENOME_END_TAG;

	/* static methods */
	/** get the default GFF annotation filename from a database name */
	static string getDBAnnoFn(const string& dbName) {
		return dbName + GFF::GFF3_SUFFIX;
	}

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

	/** read given genome info from comment */
	static istream& readStartComment(istream& in, Genome& genome);

	/** write start comment of given genome anno to text GFF output */
	static ostream& writeStartComment(ostream& out, const Genome& genome);

	/** write end comment of this genome anno to text GFF output */
	static ostream& writeEndComment(ostream& out);
};

inline MetaGenomeAnno::MetaGenomeAnno(const MetaGenome& mtg) {
	for(const Genome& genome : mtg.getGenomes())
		addGenome(genome);
}

inline bool operator==(const MetaGenomeAnno& lhs, const MetaGenomeAnno& rhs) {
	return lhs.genomeAnnos == rhs.genomeAnnos && lhs.chromAnnos == rhs.chromAnnos;
}

inline bool operator!=(const MetaGenomeAnno& lhs, const MetaGenomeAnno& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOMEANNO_H_ */
