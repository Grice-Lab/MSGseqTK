/*
 * Loc.h
 *
 *  Created on: Jun 6, 2018
 *      Author: zhengqi
 */

#ifndef SRC_LOC_H_
#define SRC_LOC_H_

#include <cassert>
#include <cstdint>
#include <iostream>

namespace EGriceLab {
namespace MSGseqTK {

using std::istream;
using std::ostream;

/** a basic location class */
class Loc {
public:
	/** default constructor */
	Loc() = default;

	/** constructor given values */
	Loc(int64_t start, int64_t end) : start(start), end(end)
	{  }

	/** destructor */
	virtual ~Loc() {  }

	/* getters and setters */
	int64_t getEnd() const {
		return end;
	}

	void setEnd(int64_t end = 0) {
		this->end = end;
	}

	int64_t getStart() const {
		return start;
	}

	void setStart(int64_t start = 0) {
		this->start = start;
	}

	/* member methods */
	int64_t length() const {
		return end - start;
	}

	bool empty() const {
		return length() <= 0;
	}

	/** reverse a Loc with given size */
	Loc& reverse(int64_t size) {
		reverseLoc(size, start, end);
		return *this;
	}

	/** get a copy of reversed GenomeLoc */
	Loc reverse(int64_t size) const {
		Loc rLoc(*this);
		return rLoc.reverse(size);
	}

	/** save this loc to binary output */
	virtual ostream& save(ostream& out) const;

	/** load a loc from a binary input */
	virtual istream& load(istream& in);

	/** write this loc to text output */
	virtual ostream& write(ostream& out) const;

	/** read this loc from text output */
	virtual istream& read(istream& in);

	/** static member methods */
	static bool isOverlap(const Loc& lhs, const Loc& rhs) {
		return lhs.start < rhs.end && lhs.end > rhs.start;
	}

	static int64_t dist(const Loc& lhs, const Loc& rhs) {
		return isOverlap(lhs, rhs) ? 0 : lhs.start < rhs.start ? rhs.start - lhs.end + 1: lhs.start - rhs.end + 1;
	}

	/**
	 * reverse a pos given total size
	 * @param size  size of the region
	 * @param i  0 or 1-based pos
	 * @return  1 or 0-based reversed pos
	 */
	static int64_t reverseLoc(int64_t size, int64_t i) {
		assert(i <= size);
		return size - i;
	}

	/** reverse a loc region with given size */
	static void reverseLoc(int64_t size, int64_t& start, int64_t& end);

	/* non-member methods */
	/** formatted output */
	friend ostream& operator<<(ostream& out, const Loc& loc);
	/** formatted input */
	friend istream& operator>>(istream& in, Loc& loc);

	/* relational operators */
	friend bool operator==(const Loc& lhs, const Loc& rhs);
	friend bool operator<(const Loc& lhs, const Loc& rhs);

private:
	/* member fields */
	int64_t start = 0; /* 0-based */
	int64_t end = 0;   /* 1-based */
};

inline void Loc::reverseLoc(int64_t size, int64_t& start, int64_t& end) {
	int64_t tmp = start;
	start = reverseLoc(size, end);
	end = reverseLoc(size, tmp);
}

inline ostream& operator<<(ostream& out, const Loc& loc) {
	return loc.write(out); /* call virtual member method */
}

inline istream& operator>>(istream& in, Loc& loc) {
	return loc.read(in);
}

inline bool operator==(const Loc& lhs, const Loc& rhs) {
	return lhs.start == rhs.start && lhs.end == rhs.end;
}

inline bool operator<(const Loc& lhs, const Loc& rhs) {
	return lhs.start != rhs.start ? lhs.start < rhs.start : lhs.end < rhs.end;
}

inline bool operator!=(const Loc& lhs, const Loc& rhs) {
	return !(lhs == rhs);
}

inline bool operator<=(const Loc& lhs, const Loc& rhs) {
	return lhs < rhs || lhs == rhs;
}

inline bool operator>(const Loc& lhs, const Loc& rhs) {
	return !(lhs <= rhs);
}

inline bool operator>=(const Loc& lhs, const Loc& rhs) {
	return !(lhs < rhs);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */

/** template specialization for customized hash function in std namespace */
namespace std {
template<>
class hash<EGriceLab::MSGseqTK::Loc> {
public:
  size_t operator() (const EGriceLab::MSGseqTK::Loc& loc) const {
	size_t res = 0;
	res ^= loc.getStart() + 0x9e3779b9 + (res << 6) + (res >> 2);
	res ^= loc.getEnd() + 0x9e3779b9 + (res << 6) + (res >> 2);
	return res;
  }
};

} /* namespace std */

#endif /* SRC_LOC_H_ */
