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
	/** default constructor */
	Loc() = default;

	/** constructor from given values */
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
	static bool isOverlap(const Loc& loc1, const Loc& loc2) {
		return loc1.start < loc2.end && loc1.end > loc2.start;
	}

	static int64_t dist(const Loc& loc1, const Loc& loc2);

	/* non-member methods */
	friend ostream& operator<<(ostream& out, const Loc& loc);
	friend istream& operator>>(istream& in, Loc& loc);

	/* member fields */
	int64_t start = 0; /* 0-based */
	int64_t end = 0;   /* 1-based */
};

inline ostream& operator<<(ostream& out, const Loc& loc) {
	return loc.write(out); /* call virtual member method */
}

inline istream& operator>>(istream& in, Loc& loc) {
	return loc.read(in);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

#endif /* SRC_LOC_H_ */
