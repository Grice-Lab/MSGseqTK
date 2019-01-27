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
	typedef vector<string> NAME_INDEX;    /* id->name vector */
	typedef map<string, size_t> GENOME_INDEX; /* genome name->id map */
	typedef map<string, size_t> CHROM_INDEX;  /* chrom name->id map */
	typedef vector<size_t> INDEX_MAP;  /* chrom idx->genome idx */
	typedef vector<size_t> BEFORE_MAP; /* chrom idx-># chroms before */
	typedef vector<Loc> GENOME_LOC; /* genome id->Loc map */
	typedef vector<Loc> CHROM_LOC;  /* chromosome id->Loc map */

	/* constructors */

	/** get total size of this MetaGenome */
	uint64_t size() const;

	/** test whether this MetaGnome is empty */
	bool empty() const {
		return genomes.empty();
	}

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

	/** get all genomeIds */
	const NAME_INDEX& getGenomeIds() const {
		return genomeIds;
	}

	/** get all chromNames */
	const NAME_INDEX& getChromNames() const {
		return chromNames;
	}

	/** get genomeId2Idx */
	const GENOME_INDEX& getGenomeId2Idx() const {
		return genomeId2Idx;
	}

	/** get chromName2Idx */
	const CHROM_INDEX& getChromName2Idx() const {
		return chromName2Idx;
	}

	/** get genomeIdx2Loc */
	const GENOME_LOC& getGenomeIdx2Loc() const {
		return genomeIdx2Loc;
	}

	/** get chromName2Loc */
	const CHROM_LOC& getChromIdx2Loc() const {
		return chromIdx2Loc;
	}

	/** get genomeId of given pos */
	const string& getGenomeId(size_t i) const {
		return genomeIds[i];
	}

	/** get genomeName of given pos */
	const string& getGenomeName(size_t i) const {
		return getGenome(i).getName();
	}

	/** get chromName of given pos */
	const string& getChromName(size_t i) const {
		return chromNames[i];
	}

	/**
	 * get genome index by id
	 * @return  index of given genome id,
	 * or -1 if not found
	 */
	size_t getGenomeIndex(const string& genomeId) const {
		return genomeId2Idx.count(genomeId) > 0 ? genomeId2Idx.at(genomeId) : -1;
	}

	/**
	 * get the genome index at given location,
	 * or -1 if not found
	 */
	size_t getGenomeIndex(uint64_t loc) const;

	/** get genome index by chrom index */
	size_t getGenomeIndexByChromIdx(size_t chrIdx) const {
		return chromIdx2GenomeIdx[chrIdx];
	}

	/**
	 * get chrom index by name
	 * @return  index of given chrom,
	 * or -1 if not found
	 */
	size_t getChromIndex(const string& chrName) const {
		return chromName2Idx.count(chrName) > 0 ? chromName2Idx.at(chrName) : -1;
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

	/** get # chroms before given chrom index in its genome */
	size_t getChromNbefore(size_t chromIdx) const {
		return chromIdx2Nbefore[chromIdx];
	}

	/** get relative chromId of a chrom index, alias as getChromNbefore */
	size_t getRelChromIndex(size_t chromIdx) const {
		return getChromNbefore(chromIdx);
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
		return genomeIdx2Loc[i];
	}

	/**
	 * get the start of the i-th genome
	 */
	int64_t getGenomeStart(size_t i) const {
		return getGenomeLoc(i).start;
	}

	/**
	 * get the end of the i-th genome, including GAP_BASE-terminal
	 */
	int64_t getGenomeEnd(size_t i) const {
		return getGenomeLoc(i).end;
	}

	/** get the length of the i-th genome, including GAP_BASE-terminal */
	uint64_t getGenomeLen(size_t i) const {
		return getGenomeLoc(i).length();
	}

	/**
	 * get the Loc of i-th chrom
	 */
	const Loc& getChromLoc(size_t i) const {
		return chromIdx2Loc[i];
	}

	/**
	 * get the start of the i-th chrom
	 */
	int64_t getChromStart(size_t i) const {
		return getChromLoc(i).start;
	}

	/**
	 * get the end of the i-th chrom, including GAP_BASE-terminal
	 */
	int64_t getChromEnd(size_t i) const {
		return getChromLoc(i).end;
	}

	/**
	 * get the length of the i-th chrom, including GAP_BASE-terminal
	 */
	uint64_t getChromLen(size_t i) const {
		return getChromLoc(i).length();
	}

	/**
	 * get the idx->Loc map for all genomes
	 */
	const GENOME_LOC& getGenomeLocs() const {
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
	void addGenome(const Genome& genome, const DNAseq& genomeSeq);

	/** get chrom seq by genome index and chrom index, including GAP_BASE terminal */
	DNAseq getChromSeq(size_t chrIdx, bool terminalGap = true) const {
		return terminalGap ? seq.substr(getChromStart(chrIdx), getChromLen(chrIdx)) :
				seq.substr(getChromStart(chrIdx), getChromLen(chrIdx) - 1);
	}

	/** get genome seq by genome index, including GAP_BASE terminal */
	DNAseq getGenomeSeq(size_t genomeIdx, bool terminalGap = true) const {
		return terminalGap ? seq.substr(getGenomeStart(genomeIdx), getGenomeLen(genomeIdx)) :
				seq.substr(getGenomeStart(genomeIdx), getGenomeLen(genomeIdx) - 1);
	}

	/** get the entire metagenome seq */
	const DNAseq& getSeq() const {
		return seq;
	}

	/** get a segment/subseq of this metagenome */
	DNAseq subseq(size_t pos = 0, size_t len = DNAseq::npos) const {
		return seq.substr(pos, len);
	}

	/** save this object to binary output, in compressed or raw mode */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in, bool ignoreSeq = false);

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
	friend MetaGenome operator+(const MetaGenome& lhs, const MetaGenome& rhs);

	/** relational operators */
	friend bool operator==(const MetaGenome& lhs, const MetaGenome& rhs);

	/* member fields */
private:
	vector<Genome> genomes;
	DNAseq seq; // GAP_BASE seperated contatenated sequence of this metagenome
	NAME_INDEX genomeIds;
	NAME_INDEX chromNames;
	GENOME_INDEX genomeId2Idx;
	CHROM_INDEX chromName2Idx;
	INDEX_MAP chromIdx2GenomeIdx;
	BEFORE_MAP chromIdx2Nbefore;
	GENOME_LOC genomeIdx2Loc; // index->0-based start
	CHROM_LOC chromIdx2Loc; // index->0-based start

public:
	/* static methods */
	/** get a unique chrom id from genome.id and chrom.name */
	static string getChromId(const string& genomeId, const string& chrName) {
		return genomeId + ":" + chrName;
	}
};

inline bool operator==(const MetaGenome& lhs, const MetaGenome& rhs) {
	return lhs.genomes == rhs.genomes && lhs.seq == rhs.seq; /* equal genomes guarantees equal index */
}

inline bool operator!=(const MetaGenome& lhs, const MetaGenome& rhs) {
	return !(lhs == rhs);
}

inline size_t MetaGenome::getGenomeIndex(uint64_t loc) const {
	GENOME_LOC::const_iterator result = std::find_if(genomeIdx2Loc.begin(), genomeIdx2Loc.end(),
			[=] (const GENOME_LOC::value_type& item) { return item.start <= loc && loc < item.end; }
	);
	return result != genomeIdx2Loc.end() ? result - genomeIdx2Loc.begin() : -1;
}

inline size_t MetaGenome::getChromIndex(uint64_t loc) const {
	CHROM_LOC::const_iterator result = std::find_if(chromIdx2Loc.begin(), chromIdx2Loc.end(),
			[=] (const CHROM_LOC::value_type& item) { return item.start <= loc && loc < item.end; }
	);
	return result != chromIdx2Loc.end() ? result - chromIdx2Loc.begin() : -1;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
