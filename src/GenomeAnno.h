/*
 * GenomeAnno.h
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */

#ifndef SRC_GENOMEANNO_H_
#define SRC_GENOMEANNO_H_

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "Genome.h"
#include "GFF.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::string;
using std::vector;
using std::map;
using std::istream;
using std::ostream;
using UCSC::GFF;

/**
 *A annotation class for a given genome that stores external or internal GFF annotations
 */
class GenomeAnno {
public:
	typedef map<string, vector<GFF>> CHROM_ANNOMAP;

	/* constructors */
	/** default constructor */
	GenomeAnno() = default;

	/** construct a GenomeAnno with empty annotations */
	GenomeAnno(const Genome& genome) : genome(genome)
	{  }

	/* member methods */
	/** get number of annotated chromosomes */
	size_t numAnnotated() const {
		return chromAnnos.size();
	}

	/** add a new GFF record for given genomeID */
	GenomeAnno& addAnno(const string& chr, const GFF& record) {
		chromAnnos[chr].push_back(record);
		return *this;
	}

	/** save to binary output */
	ostream& save(ostream& out) const;

	/** load from binary input */
	istream& load(istream& in);

	/** write MetaGenome annotation to GFF3 file in internal order */
	ostream& write(ostream& out) const;

	/** read into this object from text GFF format */
	istream& read(istream& in, GFF::Version ver = FORMAT);

	/* member fields */
	Genome genome; /* genome these annotations belong to */
	CHROM_ANNOMAP chromAnnos;

	/* static member fields */
	static const GFF::Version FORMAT = GFF::GFF3;

	/* static methods */
	/** write given genome info to comment */
	static ostream& writeGFFComment(ostream& out, const string& id, const string& name);

	/** read given genome info from comment */
	static istream& readGFFComment(istream& in, string& id, string& name);
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_GENOMEANNO_H_ */
