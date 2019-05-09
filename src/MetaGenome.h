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
#include <cassert>
#include "DNAseq.h"
#include "Genome.h"
#include "GLoc.h"
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
 * class representing a MetaGenome,
 * with its basic information and sequences stored in separately files
 */
class MetaGenome {
public:
	typedef vector<string> NAME_INDEX;    /* id->name vector */
	typedef map<string, size_t> GENOME_INDEX; /* genome id->idx map */
	typedef map<string, size_t> CHROM_INDEX;  /* chrom name->idx map */
	typedef vector<size_t> INDEX_MAP;  /* idx->id map */
	typedef vector<Loc> CHROM_LOC;  /* chromosome id->Loc map */

	/** get size of this MetaGenome, include GAP_BASE flanking each chrom */
	size_t size() const;

	/** get bi-directional size of this MetaGenome, including GAP_BASE flanking every fwd/rev chrom */
	size_t BDSize() const {
		return size() * 2;
	}

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
//	const GENOME_LOC& getGenomeIdx2Loc() const {
//		return genomeIdx2Loc;
//	}

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
	size_t getGenomeIndex(const string& genomeId) const;

	/**
	 * get the genome index at given position,
	 * or -1 if not found
	 */
	size_t getGenomeIndex(int64_t pos) const;

	/** get genome index by chrom index */
	size_t getGenomeIndexByChromIdx(size_t chrIdx) const {
		return chromIdx2GenomeIdx[chrIdx];
	}

	/**
	 * get chrom index by name
	 * @return  index of given chrom,
	 * or -1 if not found
	 */
	size_t getChromIndex(const string& chrName) const;

	/**
	 * get tid by BD position
	 */
	int64_t getLocId(int64_t pos) const;

	/** get loc on this MetaGenome with given tid and BDLoc index */
	size_t getLoc(size_t tid, int64_t pos) const {
		const Loc& tLoc = getChromBDLoc(tid);
		assert(tLoc.getStart() <= pos && pos < tLoc.getEnd());
		if(getStrand(tid, pos) == GLoc::FWD)
			return getChromLoc(tid).getStart() + pos - tLoc.getStart(); // dist from fwd start to pos
		else
			return getChromLoc(tid).getStart() + tLoc.getEnd() - pos - 1; // dist from rev end to pos
	}

	/** get loc on this Metagenome given BD pos*/
	size_t getLoc(int64_t pos) const {
		return getLoc(getLocId(pos), pos);
	}

	/** get strand of given chrom index and position seq */
	GLoc::STRAND getStrand(size_t tid, int64_t pos) const {
		const Loc& bdLoc = getChromBDLoc(tid);
		assert(bdLoc.getStart() <= pos && pos < bdLoc.getEnd());
		assert(bdLoc.length() % 2 == 0);
		return pos < bdLoc.getStart() + bdLoc.length() / 2 ? GLoc::FWD : GLoc::REV;
	}

	/**
	 * get the strand of given BD position
	 */
	GLoc::STRAND getStrand(int64_t pos) const {
		return getStrand(getLocId(pos), pos);
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
	 * get the genome at given position
	 */
	const Genome& getGenomeAtLoc(int64_t pos) const {
		return genomes.at(getGenomeIndex(pos));
	}

	/**
	 * get the Loc of i-th genome
	 */
//	const Loc& getGenomeLoc(size_t i) const {
//		return genomeIdx2Loc[i];
//	}

	/**
	 * get the start of the i-th genome
	 */
//	int64_t getGenomeStart(size_t i) const {
//		return getGenomeLoc(i).getStart();
//	}

	/**
	 * get the end of the i-th genome
	 */
//	int64_t getGenomeEnd(size_t i) const {
//		return getGenomeLoc(i).getEnd();
//	}

	/** get the length of the i-th genome, include fwd + rev */
//	int64_t getGenomeLen(size_t i) const {
//		return getGenomeLoc(i).length();
//	}

	/**
	 * get the Loc of i-th chrom
	 */
	const Loc& getChromLoc(size_t i) const {
		return chromIdx2Loc[i];
	}

	/**
	 * get the BD Loc of i-th chrom
	 */
	Loc getChromBDLoc(size_t i) const {
		return chromIdx2BDLoc[i];
	}

	/** get the BD Loc of a region */
	Loc getChromBDLoc(size_t tidStart, size_t tidEnd) const {
		assert(tidStart < tidEnd);
		return Loc(chromIdx2BDLoc[tidStart].getStart(), chromIdx2BDLoc[tidEnd - 1].getEnd());
	}

	/**
	 * get the idx->Loc map for all genomes
	 */
//	const GENOME_LOC& getGenomeLocs() const {
//		return genomeIdx2Loc;
//	}

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
	void addGenome(const Genome& genome) {
		genomes.push_back(genome);
	}

//	/** get the last/top genome, const version */
//	const Genome& topGenome() const {
//		return genomes.back();
//	}
//
//	/** get the last/top genome, non-const version */
//	Genome& topGenome() {
//		return genomes.back();
//	}
//
//	/** pop the last/top genome */
//	void popGenome() {
//		genomes.pop_back();
//	}
//
//	/** get the last/top chromsome, const version */
//	const Genome::Chrom& topChrom() const {
//		return genomes.back().chroms.back();
//	}
//
//	/** get the last/top chromsome, non-const version */
//	Genome::Chrom& topChrom() {
//		return genomes.back().chroms.back();
//	}
//
//	/** pop the last chorosome */
//	void popChrom() {
//		genomes.back().chroms.pop_back();
//		if(genomes.back().chroms.empty()) // no more chromsomes
//			genomes.pop_back();
//	}

	/** get chromosome given chrIdx */
	const Genome::Chrom& getChrom(size_t chrIdx) const {
		return getGenome(getGenomeIndexByChromIdx(chrIdx)).getChrom(getChromNbefore(chrIdx));
	}

	/** get seq of the entire metagenome */
	const DNAseq& getSeq() const {
		return seq;
	}

	/** get a chrom seq by tid */
	DNAseq getSeq(size_t tid) const {
		return seq.substr(getChromLoc(tid).getStart(), getChromLoc(tid).length());
	}

	/** get a block seq by tid range */
	DNAseq getSeq(size_t tidStart, size_t tidEnd) const {
		assert(tidStart < tidEnd);
		return seq.substr(getChromLoc(tidStart).getStart(), getChromLoc(tidEnd - 1).getEnd() - getChromLoc(tidStart).getStart());
	}

	/** get bi-directional seq by tid */
	DNAseq getBDSeq(size_t tid) const {
		return getBDSeq(getSeq(tid));
	}

	/** get a BD block seq by tid range */
	DNAseq getBDSeq(size_t tidStart, size_t tidEnd) const {
		assert(tidStart < tidEnd);
		DNAseq seq;
		seq.reserve(getChromBDLoc(tidEnd - 1).getEnd() - getChromBDLoc(tidStart).getStart());
		for(size_t tid = tidStart; tid < tidEnd; ++tid)
			seq += getBDSeq(tid);
		return seq;
	}

	/** get concatenated bi-directional seq of the entire MetaGenome */
	DNAseq getBDSeq() const {
		return getBDSeq(0, numChroms());
	}

	/** remove seq of tid, will invalid the tid unless remove backwards */
	void removeSeq(size_t tid) {
		seq.erase(getChromLoc(tid).getStart(), getChromLoc(tid).length());
	}

	/** remove a block */
	void removeSeq(size_t tidStart, size_t tidEnd) {
		assert(tidStart < tidEnd);
		seq.erase(getChromLoc(tidStart).getStart(), getChromLoc(tidEnd - 1).getEnd() - getChromLoc(tidStart).getStart());
	}

	/** get gap location in the BDseq */
	vector<size_t> getBDGapLoc() const;

	/** save this object to binary outputs */
	ostream& save(ostream& out) const;

	/** load an object from binary inputs */
	istream& load(istream& in);

	/** update MetaGenome */
	void update() {
		updateIndex();
		updateSeq();
	}

	/**
	 * update seq by moving seqs from Chrom into independent storage
	 */
	void updateSeq();

	/**
	 * update all index, should be called if any containing Genome/Chrom changes
	 * it will rename genome ids and chrom names (with warnings), if non-unique ID/names found
	 */
	void updateIndex();

	/** clear stored seq */
	void clearSeq() {
		seq.clear();
	}

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
	DNAseq seq; // conncatenated metagenome seq (forward strand only, each chrom separated by GAP_BASE
	NAME_INDEX genomeIds; // all genome ids
	NAME_INDEX chromNames; // all chrom names
	GENOME_INDEX genomeId2Idx; // genome id -> idx
	CHROM_INDEX chromName2Idx; // chrom name -> idx
	INDEX_MAP chromIdx2GenomeIdx; // chrom idx -> genome idx
	INDEX_MAP chromIdx2Nbefore; // chrom idx -> # chroms before this chrom in this genome
//	GENOME_LOC genomeIdx2Loc; // genome idx -> Loc on MetaGenome (including fwd + rev)
	CHROM_LOC chromIdx2Loc; // chrom idx -> Loc on MetaGenome
	CHROM_LOC chromIdx2BDLoc; // chrom idx -> bi-directional Loc on MetaGenome

public:
	/* static methods */
	/** get a unique chrom id from genome.id and chrom.name */
	static string getChromId(const string& genomeId, const string& chrName) {
		return genomeId + ":" + chrName;
	}

	/** get bi-directional seq for a given seq */
	static DNAseq getBDSeq(const DNAseq& seq) {
		return dna::toBasic(seq) + dna::revcom(dna::toBasic(seq));
	}
};

inline bool operator==(const MetaGenome& lhs, const MetaGenome& rhs) {
	return lhs.genomes == rhs.genomes && lhs.seq == rhs.seq;
}

inline bool operator!=(const MetaGenome& lhs, const MetaGenome& rhs) {
	return !(lhs == rhs);
}

inline size_t MetaGenome::getGenomeIndex(const string& genomeId) const {
	GENOME_INDEX::const_iterator result = genomeId2Idx.find(genomeId);
	return result != genomeId2Idx.end() ? result->second : -1;
}

//inline size_t MetaGenome::getGenomeIndex(int64_t pos) const {
//	GENOME_LOC::const_iterator result = std::find_if(genomeIdx2Loc.begin(), genomeIdx2Loc.end(),
//			[=] (const GENOME_LOC::value_type& item) { return item.getStart() <= pos && pos < item.getEnd(); }
//	);
//	return result != genomeIdx2Loc.end() ? result - genomeIdx2Loc.begin() : -1;
//}

inline size_t MetaGenome::getChromIndex(const string& chrName) const {
	CHROM_INDEX::const_iterator result = chromName2Idx.find(chrName);
	return result != chromName2Idx.end() ? result->second : -1;
}

inline int64_t MetaGenome::getLocId(int64_t pos) const {
	CHROM_LOC::const_iterator result = std::find_if(chromIdx2BDLoc.begin(), chromIdx2BDLoc.end(),
			[=] (const CHROM_LOC::value_type& item) { return item.getStart() <= pos && pos < item.getEnd(); }
	);
	return result != chromIdx2BDLoc.end() ? result - chromIdx2BDLoc.begin() : -1;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
