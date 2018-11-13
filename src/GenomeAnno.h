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
	/** get annotating genome */
	const Genome& getGenome() const {
		return genome;
	}

	/** set annotating genome */
	void setGenome(const Genome& genome) {
		this->genome = genome;
	}

	/** get number of annotated chromosomes */
	size_t numChromAnnotated() const {
		return chromAnnos.size();
	}

	/** get number of total annotations */
	size_t numAnnotated() const;

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

	/** write start comment of this genome anno to text GFF output */
	ostream& writeStartComment(ostream& out) const;

	/** write end comment of this genome anno to text GFF output */
	ostream& writeEndComment(ostream& out) const;

	/* member fields */
private:
	Genome genome; /* genome these annotations belong to */
	CHROM_ANNOMAP chromAnnos;

public:
	/* static member fields */
	static const GFF::Version FORMAT = GFF::GFF3;
	static const string RECORD_START_TAG;
	static const string RECORD_END_TAG;

	/* static methods */
	/** read given genome info from comment */
	static istream& readStartComment(istream& in, string& id, string& name);

};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_GENOMEANNO_H_ */
