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
namespace MSGseqClean {

using std::istream;
using std::ostream;

struct Loc {
	/* default constructor */
	Loc() = default;

	/* constructor from given values */
	Loc(uint64_t start, uint64_t end) : start(start), end(end)
	{  }

	/* member methods */
	uint64_t length() const {
		return end - start;
	}

	/** save this loc to binary output */
	ostream& save(ostream& out) const;

	/** load a loc from a binary input */
	istream& load(istream& in);

	friend std::ostream& operator<<(std::ostream& out, const Loc& loc);
	friend std::istream& operator>>(std::istream& in, Loc& loc);

	uint64_t start; /* 0-based */
	uint64_t end;   /* 1-based */
};

inline std::ostream& operator<<(std::ostream& out, const Loc& loc) {
	out << "[" << loc.start << ", " << loc.end << ")"; /* [start, end) */
	return out;
}

inline ostream& Loc::save(ostream& out) const {
	out.write((const char*) &start, sizeof(uint64_t));
	out.write((const char*) &end, sizeof(uint64_t));
	return out;
}

inline istream& Loc::load(istream& in) {
	in.read((char*) &start, sizeof(uint64_t));
	in.read((char*) &end, sizeof(uint64_t));
	return in;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */

#endif /* SRC_LOC_H_ */
