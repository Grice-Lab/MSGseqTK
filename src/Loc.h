/*
 * Loc.h
 *
 *  Created on: Jun 6, 2018
 *      Author: zhengqi
 */

#ifndef SRC_LOC_H_
#define SRC_LOC_H_

#include <cstdint>
#include <iostream>

namespace EGriceLab {
namespace MSGseqTK {

using std::istream;
using std::ostream;

struct Loc {
	/* enums */
	enum STRAND { FWD = 1 /* forward */, REV /* reverse-complement */, UNK /* unknown */  };

	/** default constructor */
	Loc() = default;

	/** constructor given all values */
	Loc(int64_t start, int64_t end, STRAND strand) : start(start), end(end), strand(strand)
	{  }

	/** constructor given coordinates */
	Loc(int64_t start, int64_t end) : start(start), end(end)
	{  }

	/* member methods */
	uint64_t length() const {
		return end - start;
	}

	/** save this loc to binary output */
	ostream& save(ostream& out) const;

	/** load a loc from a binary input */
	istream& load(istream& in);

	/** write this loc to text output */
	ostream& write(ostream& out) const;

	/** read this loc from text output */
	istream& read(istream& in);

	/** static member methods */
	static bool isOverlap(const Loc& lhs, const Loc& rhs, bool ignoreStrand = true) {
		return lhs.start < rhs.end && lhs.end > rhs.start &&
				(ignoreStrand || (lhs.strand & rhs.strand) != 0);
	}

	static int64_t dist(const Loc& lhs, const Loc& rhs) {
		return isOverlap(lhs, rhs) ? 0 : lhs.start < rhs.start ? rhs.start - lhs.end + 1: lhs.start - rhs.end + 1;
	}

	/** decode strand to char */
	static char decodeStrand(STRAND strand);

	/** encode strand from char */
	static STRAND encodeStrand(char s);

	/* non-member methods */
	/** formatted output */
	friend ostream& operator<<(ostream& out, const Loc& loc);
	/** formatted input */
	friend istream& operator>>(istream& in, Loc& loc);

	/* relational operators */
	friend bool operator==(const Loc& lhs, const Loc& rhs);

	/* member fields */
	int64_t start = 0; /* 0-based */
	int64_t end = 0;   /* 1-based */
	STRAND strand = UNK;
};

inline ostream& operator<<(ostream& out, const Loc& loc) {
	return loc.write(out); /* call virtual member method */
}

inline istream& operator>>(istream& in, Loc& loc) {
	return loc.read(in);
}

inline bool operator==(const Loc& lhs, const Loc& rhs) {
	return lhs.start == rhs.start && lhs.end == rhs.end;
}

inline bool operator!=(const Loc& lhs, const Loc& rhs) {
	return !(lhs == rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_LOC_H_ */
