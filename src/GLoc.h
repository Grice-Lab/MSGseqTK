/*
 * GLoc.h
 *
 *  Created on: Feb 13, 2019
 *      Author: zhengqi
 */

#ifndef SRC_GLOC_H_
#define SRC_GLOC_H_

#include "Loc.h"

namespace EGriceLab {
namespace MSGseqTK {

/** an extended location class with chrom and strand info */
class GLoc: public Loc {
public:
	/* nested enums */
	enum STRAND { FWD = 1 /* forward */, REV /* reverse-complement */, UNK /* unknown */  };

	/* constructors */
	/** default constructor */
	GLoc() = default;

	/** construct with all given values */
	GLoc(int64_t start, int64_t end, int64_t tid, STRAND strand) : Loc(start, end), tid(tid), strand(strand)
	{  }

	/** construct with basic values */
	GLoc(int64_t start, int64_t end) : Loc(start, end)
	{  }

	virtual ~GLoc() {  }

	/* getters and setters */
	STRAND getStrand() const {
		return strand;
	}

	void setStrand(STRAND strand = UNK) {
		this->strand = strand;
	}

	int64_t getTid() const {
		return tid;
	}

	void setTid(int64_t tid = -1) {
		this->tid = tid;
	}

	/* member methods */
	/** complement a GenomeLoc */
	GLoc& complement() {
		switch(strand) {
		case FWD:
			strand = REV;
			break;
		case REV:
			strand = FWD;
			break;
		}
		return *this;
	}

	/** get a copy of complement GenomeLoc */
	GLoc complement() const {
		GLoc cLoc(*this);
		return cLoc.complement();
	}

	/** reverse-complement a GenomeLoc with given size */
	GLoc& revcom(int64_t size) {
		reverse(size);
		return complement();
	}

	/** get a copy of reverse-complement loc */
	GLoc revcom(int64_t size) const {
		GLoc rcLoc(*this);
		return rcLoc.revcom(size);
	}

	/**
	 * save this loc to binary output
	 * @override base class method
	 */
	virtual ostream& save(ostream& out) const;

	/**
	 * load a loc from a binary input
	 * @override base class method
	 */
	virtual istream& load(istream& in);

	/**
	 * write this loc to text output
	 * @override base class method
	 */
	virtual ostream& write(ostream& out) const;

	/**
	 * read this loc from text output
	 * @override base class method
	 */
	virtual istream& read(istream& in);

	/** static member methods */
	static bool isOverlap(const GLoc& lhs, const GLoc& rhs) {
		return lhs.tid == rhs.tid && lhs.strand == rhs.strand && Loc::isOverlap(lhs, rhs);
	}

	static int64_t dist(const GLoc& lhs, const GLoc& rhs) {
		return isOverlap(lhs, rhs) ? 0 : Loc::dist(lhs, rhs);
	}

	/** decode strand to char */
	static char decodeStrand(STRAND strand);

	/** encode strand from char */
	static STRAND encodeStrand(char s);

	/* non-member methods */
	/* relational operators */
	friend bool operator==(const GLoc& lhs, const GLoc& rhs);
	friend bool operator<(const GLoc& lhs, const GLoc& rhs);

	/* member fields */
private:
	int64_t tid = -1;
	STRAND strand = UNK;
};

inline bool operator==(const GLoc& lhs, const GLoc& rhs) {
	return dynamic_cast<const Loc&>(lhs) == dynamic_cast<const Loc&>(rhs)
			&& lhs.tid == rhs.tid && lhs.strand == rhs.strand;
}

inline bool operator<(const GLoc& lhs, const GLoc& rhs) {
	return dynamic_cast<const Loc&>(lhs) != dynamic_cast<const Loc&>(rhs) ?
			dynamic_cast<const Loc&>(lhs) < dynamic_cast<const Loc&>(rhs) :
			lhs.tid != rhs.tid ? lhs.tid < rhs.tid : lhs.strand < rhs.strand;
}

inline bool operator!=(const GLoc& lhs, const GLoc& rhs) {
	return !(lhs == rhs);
}

inline bool operator<=(const GLoc& lhs, const GLoc& rhs) {
	return lhs < rhs || lhs == rhs;
}

inline bool operator>(const GLoc& lhs, const GLoc& rhs) {
	return rhs < lhs;
}

inline bool operator>=(const GLoc& lhs, const GLoc& rhs) {
	return !(lhs < rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

/** template specialization for customized hash function in std namespace */
namespace std {
template<>
class hash<EGriceLab::MSGseqTK::GLoc> {
public:
  size_t operator() (const EGriceLab::MSGseqTK::GLoc& gLoc) const {
	size_t res = 0;
	res ^= gLoc.getStart() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= gLoc.getEnd() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= gLoc.getTid() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= gLoc.getStrand() + 0x9e3779b9 + (res << 6) + (res >> 2);
	return res;
  }
};

} /* namespace std */

#endif /* SRC_GLOC_H_ */
