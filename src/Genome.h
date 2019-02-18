/*
 * Genome.h
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#ifndef SRC_GENOME_H_
#define SRC_GENOME_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cstdint> // C++11
#include <regex>
#include "DNAseq.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::map;
using std::vector;
using std::istream;
using std::ostream;

/**
 * Basic information of a genome (within a MetaGenome)
 */
class Genome {
public:
	/* nested types and enums */
	struct Chrom {
		/** default constructor */
		Chrom() = default;

		/** construct from given info */
		Chrom(const string& name, const DNAseq& seq) : name(name), seq(seq)
		{  }

		/* member methods */
		size_t size() const {
			return seq.length();
		}

		/** save this Chrom to binary output */
		ostream& save(ostream& out) const;

		/** load a Chrom from binary input */
		istream& load(istream& in);

		/* non-member functions */
		friend bool operator==(const Chrom& lhs, const Chrom& rhs);

		string name;
		DNAseq seq;
	};

	/* constructors */
	/** default constructor */
	Genome() = default;

	/** construct genome with given name */
	Genome(const string& id, const string& name) : id(id), name(name)
	{  }

	/** deligating construct from name only */
	Genome(const string& name) : Genome(name, name)
	{  }

	/** construct a genome with given id, name and chromosomes */
	Genome(const string& id, const string& name, const vector<Chrom>& chroms) : id(id), name(name), chroms(chroms)
	{  }

	/** deligating construct a genome with give name and chromosomes */
	Genome(const string& name, const vector<Chrom>& chroms) : Genome(name, name, chroms)
	{  }

	/* member methods */
	const string& getId() const {
		return id;
	}

	void setId(const string& id) {
		this->id = id;
	}

	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	/** get a display id of this genome */
	string displayId() const {
		return displayId(id, name);
	}

	/** get number of chromosomes */
	size_t numChroms() const {
		return chroms.size();
	}

	/** get all chromosomes */
	const vector<Chrom>& getChroms() const {
		return chroms;
	}

	/** get chromosome by index */
	const Chrom& getChrom(size_t i) const {
		return chroms[i];
	}

	/** test whether this chrom name exists */
	bool hasChrom(const string& chrName) const;

	/** get a single chromsome size, or 0 if not found */
	size_t getChromSize(const string& chrName) const;

	/** get the overall size of this genome */
	size_t size() const;

	/** get bi-directional seq */
	DNAseq getSeq() const;

	/** add a new chrom object at the end */
	void addChrom(const Chrom& chr) {
		chroms.push_back(chr);
	}

	/** add a new chromosome with given name and size */
	void addChrom(const string& chrName, const DNAseq& chrSeq) {
		addChrom(Chrom(chrName, chrSeq));
	}

	/**
	 * get chromosome index given a relative loc of this Genome
	 * @param  i  0-based position relative to this Genome
	 * @return  0-based relative order, or -1 if not exists
	 */
	size_t getChromIndex(int64_t i) const;

	/** save this object to binary output */
	ostream& save(ostream& out) const;

	/** load an object from binary input */
	istream& load(istream& in);

	/* non-member functions */
	/* relationship operators */
	/** test whether two Genomes are equal, all basic information but not auxinary annotations are tested */
	friend bool operator==(const Genome& lhs, const Genome& rhs);

	/* member fields */
private:
	string id; /* unique genome ID useful for upgrade */
	string name;
	vector<Chrom> chroms;

	/* class members */
public:
	static const std::regex INVALID_NAMEPREFIX_PATTERN;
	static const std::regex INVALID_NAME_PATTERN;
	static const string REPLACEMENT_STR;

	/* static methods */
	/** format genome id/name */
	static string formatName(const string& name);

	/** get display id of given id and name */
	static string displayId(const string& id, const string& name) {
		return id + " (" + name + ")";
	}

	/* friend declarations */
	friend class MetaGenome;
	friend class GenomeAnno;
	friend class MetaGenomeAnno;
};

inline bool operator==(const Genome& lhs, const Genome& rhs) {
	return lhs.id == rhs.id && lhs.name == rhs.name && lhs.chroms == rhs.chroms;
}

inline bool operator==(const Genome::Chrom& lhs, const Genome::Chrom& rhs) {
	return lhs.name == rhs.name && lhs.seq == rhs.seq;
}

inline bool operator!=(const Genome::Chrom& lhs, const Genome::Chrom& rhs) {
	return !(lhs == rhs);
}

inline bool operator!=(const Genome& lhs, const Genome& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_GENOME_H_ */
