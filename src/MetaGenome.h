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
#include <algorithm>
#include "DNAseq.h"
#include "Genome.h"
#include "Loc.h"
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
	typedef map<string, size_t> GENOME_INDEX; /* genome name->id map */
	typedef map<string, size_t> CHROM_INDEX;  /* chrom name->id map */
	typedef map<size_t, Loc> GENOME_LOCMAP; /* genome id->Loc map */
	typedef map<size_t, Loc> CHROM_LOCMAP;  /* chromosome id->Loc map */

	/* constructors */

	/** get total size of this MetaGenome */
	uint64_t size() const;

	/** get total number of genomes */
	size_t numGenomes() const {
		return genomes.size();
	}

	/** get total number of chromosomes */
	size_t numChroms() const;

	/** get all Genomes in this MetaGenome */
	const vector<Genome>& getGenomes() const {
		return genomes;
	}

	/**
	 * get genome index by name
	 * @return  index of given genome name,
	 * or -1 if not found
	 */
	size_t getGenomeIndex(const string& gname) const {
		return genomeId2Idx.count(gname) > 0 ? genomeId2Idx.at(gname) : -1;
	}

	/**
	 * get the genome index at given location,
	 * or -1 if not found
	 */
	size_t getGenomeIndex(uint64_t loc) const;

	/**
	 * get chrom index by name
	 * @return  index of given chrom,
	 * or -1 if not found
	 */
	size_t getChromIndex(const string& cname) const {
		return chromName2Idx.count(cname) > 0 ? chromName2Idx.at(cname) : -1;
	}

	/**
	 * get the chromosome index of given location,
	 * or -1 if not found
	 */
	size_t getChromIndex(uint64_t loc) const;

	/**
	 * get a general id for given loc,
	 * alias to getChromIndex
	 */
	size_t getLocId(uint64_t loc) const {
		return getChromIndex(loc);
	}

	/**
	 * get the genome at given location
	 */
	const Genome& getGenomeAtLoc(uint64_t loc) const {
		return genomes.at(getGenomeIndex(loc));
	}

	/**
	 * get the Loc of i-th genome
	 */
	const Loc& getGenomeLoc(size_t i) const {
		return genomeIdx2Loc.at(i);
	}

	/**
	 * get the start of the i-th genome
	 */
	int64_t getGenomeStart(size_t i) const {
		return getGenomeLoc(i).start;
	}

	/**
	 * get the end of the i-th genome
	 */
	int64_t getGenomeEnd(size_t i) const {
		return getGenomeLoc(i).end;
	}

	/**
	 * get the Loc of i-th chrom
	 */
	const Loc& getChromLoc(size_t i) const {
		return chromIdx2Loc.at(i);
	}

	/**
	 * get the start of the i-th chrom
	 */
	int64_t getChromStart(size_t i) const {
		return getChromLoc(i).start;
	}

	/**
	 * get the end of the i-th chrom
	 */
	int64_t getChromEnd(size_t i) const {
		return getChromLoc(i).end;
	}

	/**
	 * get the idx->Loc map for all genomes
	 */
	const GENOME_LOCMAP& getGenomeLocs() const {
		return genomeIdx2Loc;
	}

	/** get genome by index */
	const Genome& getGenome(size_t i) const {
		return genomes[i];
	}

	/** get genome by name */
	const Genome& getGenome(const string& gname) const {
		return getGenome(getGenomeIndex(gname));
	}

	/** get genome by index, non-const version */
	Genome& getGenome(size_t i) {
		return genomes[i];
	}

	/** check whether this genome with given ID exists */
	bool hasGenome(const string& genomeId) const {
		return genomeId2Idx.count(genomeId) > 0;
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
		genomes.insert(genomes.begin(), genome);
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

	/**
	 * update all index, should be called if any containing Genome/Chrom changes
	 * it will rename genome ids and chrom names (with warnings), if non-unique ID/names found
	 */
	void updateIndex();

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
	vector<Genome> genomes;
	GENOME_INDEX genomeId2Idx;
	CHROM_INDEX chromName2Idx;
	GENOME_LOCMAP genomeIdx2Loc; // index->0-based start
	CHROM_LOCMAP chromIdx2Loc; // index->0-based start

public:
	/* static methods */
	/** get a unique chrom id from genome.id and chrom.name */
	static string getChromId(const string& genomeId, const string& chrName) {
		return genomeId + "." + chrName;
	}
};

inline bool operator==(const MetaGenome& lhs, const MetaGenome& rhs) {
	return lhs.genomes == rhs.genomes; /* equal genomes guarantees equal index */
}

inline bool operator!=(const MetaGenome& lhs, const MetaGenome& rhs) {
	return !(lhs == rhs);
}

inline size_t MetaGenome::getGenomeIndex(uint64_t loc) const {
	GENOME_LOCMAP::const_iterator result = std::find_if(genomeIdx2Loc.begin(), genomeIdx2Loc.end(),
			[=] (const GENOME_LOCMAP::value_type& item) { return item.second.start <= loc && loc < item.second.end; }
	);
	return result != genomeIdx2Loc.end() ? result->first : -1;
}

inline size_t MetaGenome::getChromIndex(uint64_t loc) const {
	CHROM_LOCMAP::const_iterator result = std::find_if(chromIdx2Loc.begin(), chromIdx2Loc.end(),
			[=] (const CHROM_LOCMAP::value_type& item) { return item.second.start <= loc && loc < item.second.end; }
	);
	return result != chromIdx2Loc.end() ? result->first : -1;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
