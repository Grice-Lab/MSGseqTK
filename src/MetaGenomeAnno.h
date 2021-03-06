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
#include "ProgEnv.h"
#include "GFF.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::istream;
using std::ostream;
using std::vector;
using UCSC::GFF;

/**
 * class providing static methods for MetaGenome annotation
 */
class MetaGenomeAnno {
public:
	/* static member fields */
	static const GFF::Version ANNO_VER = GFF::GFF3;
	static const string GENOME_START_TAG;
	static const string GENOME_END_TAG;
	static const string METAGENOME_ID_TAG;
	static const string METAGENOME_NAME_TAG;
	static const string METAGENOME_TYPE_TAG;
	static const string EXTERNAL_ID_TAG;
	static const string EXTERNAL_NAME_TAG;

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
	 * write metagenome GFF header
	 */
	static ostream& writeGFFHeader(ostream& out, const string& dbName, GFF::Version ver = GFF::GFF3);

	/**
	 * read pre-built GFF header comments
	 */
	static istream& readGFFHeader(istream& in, string& dbName, GFF::Version& ver);

	/** read genome from comment */
	static istream& readStartComment(istream& in, Genome& genome);

	/** write genome to GFF comment */
	static ostream& writeStartComment(ostream& out, const Genome& genome);

	/** write genome end comment */
	static ostream& writeEndComment(ostream& out);

	/** get a genome GFF annotation */
	static GFF getAnno(const Genome& genome) {
		return GFF(genome.id, progName, "genome", 1, genome.size(),
				GFF::INVALID_SCORE, GFF::DEFAULT_STRAND, GFF::INVALID_FRAME,
				GFF::attr_map { { "ID", genome.id }, { "Name", genome.name }});
	}

	/** get a chrom GFF annotation record */
	static GFF getAnno(const Genome& genome, const Genome::Chrom& chr) {
		return GFF(chr.name, progName, "chromosome", 1, chr.size(),
				GFF::INVALID_SCORE, GFF::DEFAULT_STRAND, GFF::INVALID_FRAME,
				GFF::attr_map {
			{ "ID", MetaGenome::getChromId(genome.id, chr.name) },
			{ "Name", chr.name },
			{ "Parent", genome.id },
			{ "genomeName", genome.name }
		});
	}

	/** read in all GFF records from text input in given GFF version */
	static vector<GFF> read(istream& in, GFF::Version ver);

	/**
	 * add aditional metagenome tags such as ID, Name and Type to all GFF records
	 * record with Parent feature will be traced back recursively, assuming the records are in correct order
	 */
	static vector<GFF>& addMetagenomeTags(const Genome& genome, vector<GFF>& gffRecords,
			const string& idTag = EXTERNAL_ID_TAG, const string& nameTag = EXTERNAL_NAME_TAG);

	/**
	 * write genome annotations including both genome-level and any auxilary GFF records
	 * @return  # of GFF records written
	 */
	static size_t writeGenomeAnnos(ostream& out, const Genome& genome, const vector<GFF>& gffRecords);
};

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOMEANNO_H_ */
