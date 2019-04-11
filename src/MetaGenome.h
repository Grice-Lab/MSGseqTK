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
	typedef vector<size_t> INDEX_MAP;  /* chrom idx->genome idx */
	typedef vector<size_t> BEFORE_MAP; /* chrom idx-># chroms before */
	typedef vector<size_t> STREAMPOS_MAP; /* chrom idx->stream end-pos */
//	typedef vector<Loc> GENOME_LOC; /* genome id->Loc map */
	typedef vector<Loc> CHROM_LOC;  /* chromosome id->Loc map */

	/* constructors */

	/** get size of this MetaGenome */
	size_t size() const;

	/** get bi-directional size of this MetaGenome */
	size_t BDSize() const {
		return (size() + numChroms()) * 2; // MetaGenome includes a BAP_BASE after each chromosome
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
	 * get the chromosome index of given position,
	 * or -1 if not found
	 */
	size_t getChromIndex(int64_t pos) const;

	/**
	 * get a general id for given position,
	 * alias to getChromIndex
	 */
	size_t getLocId(int64_t pos) const {
		return getChromIndex(pos);
	}

	/** get strand of given chrom index and position */
	GLoc::STRAND getStrand(size_t tid, int64_t pos) const {
		assert(getChromFwdStart(tid) <= pos && pos < getChromEnd(tid));
		return pos < getChromFwdEnd(tid) ? GLoc::FWD : GLoc::REV;
	}

	/**
	 * get the strand of given position
	 */
	GLoc::STRAND getStrand(int64_t pos) const {
		assert(getChromIndex(pos) != -1);
		return getStrand(getChromIndex(pos), pos);
	}

	/** get # chroms before given chrom index in its genome */
	size_t getChromNbefore(size_t chromIdx) const {
		return chromIdx2Nbefore[chromIdx];
	}

	/** get relative chromId of a chrom index, alias as getChromNbefore */
	size_t getRelChromIndex(size_t chromIdx) const {
		return getChromNbefore(chromIdx);
	}

	/** get end pos of a chrom index */
	size_t getEndPos(size_t chromIdx) const {
		return chromIdx2endPos[chromIdx];
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
	 * get the fwd Loc of i-th chrom
	 */
	Loc getChromFwdLoc(size_t i) const {
		return Loc(getChromFwdStart(i), getChromFwdEnd(i));
	}

	/**
	 * get the rev Loc of i-th chrom
	 */
	Loc getChromRevLoc(size_t i) const {
		return Loc(getChromRevStart(i), getChromRevEnd(i));
	}

	/**
	 * get the start of the i-th chrom
	 */
	int64_t getChromStart(size_t i) const {
		return getChromLoc(i).getStart();
	}

	/**
	 * get the fwd start of the i-th chrom
	 */
	int64_t getChromFwdStart(size_t i) const {
		return getChromStart(i);
	}

	/**
	 * get the rev start of the i-th chrom
	 */
	int64_t getChromRevStart(size_t i) const {
		return getChromStart(i) + getChromLen(i) / 2;
	}

	/**
	 * get the end of the i-th chrom
	 */
	int64_t getChromEnd(size_t i) const {
		return getChromLoc(i).getEnd();
	}

	/**
	 * get the fwd end of the i-th chrom
	 */
	int64_t getChromFwdEnd(size_t i) const {
		return getChromStart(i) + getChromLen(i) / 2;
	}

	/**
	 * get the rev end of the i-th chrom
	 */
	int64_t getChromRevEnd(size_t i) const {
		return getChromEnd(i);
	}

	/**
	 * get the length of the i-th chrom, including fwd + rev
	 */
	int64_t getChromLen(size_t i) const {
		return getChromLoc(i).length();
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

	/** get chrom seq by tid */
	const DNAseq& getSeq(size_t tid) const {
		return seqs[tid];
	}

	/** get concatenated bi-directional whole MetaGenome seq */
	DNAseq getBDSeq() const;

	/** save this object to binary outputs */
	ostream& save(ostream& out) const;

	/** load an object from binary inputs */
	istream& load(istream& in);

	/** load a chrom seq from a rw binary stream for given tid */
	DNAseq loadSeq(istream& in, size_t tid) const {
		DNAseq seq;
		in.seekg(- getEndPos(tid) /* negative value */, std::ios_base::end);
		StringUtils::loadString(seq, in);
		return seq;
	}

	/** load a chrom seq from a rw binary stream on current position */
	DNAseq loadSeq(istream& in) const {
		DNAseq seq;
		StringUtils::loadString(seq, in);
		return seq;
	}

	/**
	 * load a chrom bi-directional seq from an r/w binary input given tid
	 * the loaded seq will be convert to basic-only bases so it is bi-directional
	 */
	DNAseq loadBDSeq(istream& in, size_t tid) const {
		const DNAseq& seq = dna::toBasic(loadSeq(in, tid));
		return seq + DNAalphabet::GAP_BASE + dna::revcom(seq) + DNAalphabet::GAP_BASE;
	}

	/**
	 * load a chrom bi-directional seq from an r/w binary at current pos
	 * the loaded seq will be convert to basic-only bases so it is bi-directional
	 */
	DNAseq loadBDSeq(istream& in) const {
		const DNAseq& seq = dna::toBasic(loadSeq(in));
		return seq + DNAalphabet::GAP_BASE + dna::revcom(seq) + DNAalphabet::GAP_BASE;
	}

	/** update MetaGenome */
	void update() {
		updateSeq();
		updateIndex();
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
		seqs.clear();
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
	vector<DNAseq> seqs; // chrom seqs stored at the end of the output
	NAME_INDEX genomeIds; // all genome ids
	NAME_INDEX chromNames; // all chrom names
	GENOME_INDEX genomeId2Idx; // genome id -> idx
	CHROM_INDEX chromName2Idx; // chrom name -> idx
	INDEX_MAP chromIdx2GenomeIdx; // chrom idx -> genome idx
	STREAMPOS_MAP chromIdx2endPos; // chrom idx -> stream pos relative to the end of binary output/input
	BEFORE_MAP chromIdx2Nbefore; // chrom idx -> # chroms before this chrom in this genome
//	GENOME_LOC genomeIdx2Loc; // genome idx -> Loc on MetaGenome (including fwd + rev)
	CHROM_LOC chromIdx2Loc; // chrom idx -> Loc on MetaGenome (including fwd + rev)

public:
	/* static methods */
	/** get a unique chrom id from genome.id and chrom.name */
	static string getChromId(const string& genomeId, const string& chrName) {
		return genomeId + ":" + chrName;
	}
};

inline bool operator==(const MetaGenome& lhs, const MetaGenome& rhs) {
	return lhs.genomes == rhs.genomes; /* equal genomes guarantees equal index */
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

inline size_t MetaGenome::getChromIndex(int64_t pos) const {
	CHROM_LOC::const_iterator result = std::find_if(chromIdx2Loc.begin(), chromIdx2Loc.end(),
			[=] (const CHROM_LOC::value_type& item) { return item.getStart() <= pos && pos < item.getEnd(); }
	);
	return result != chromIdx2Loc.end() ? result - chromIdx2Loc.begin() : -1;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_METAGENOME_H_ */
