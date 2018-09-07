/*
 * Loc.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: zhengqi
 */

#include "Loc.h"

namespace EGriceLab {
namespace MSGseqTK {

ostream& Loc::save(ostream& out) const {
	out.write((const char*) &start, sizeof(int64_t));
	out.write((const char*) &end, sizeof(int64_t));
	return out;
}

istream& Loc::load(istream& in) {
	in.read((char*) &start, sizeof(int64_t));
	in.read((char*) &end, sizeof(int64_t));
	return in;
}

ostream& Loc::write(ostream& out) const {
	out << "[" << start << ", " << end << ")"; /* [start, end) */
	return out;
}

istream& Loc::read(istream& in) {
	in.ignore(1, '[');
	in >> start;
	in.ignore(1, ',');
	in >> end; /* automatical ignore additional space */
	in.ignore(1, ')');
	return in;
}

int64_t Loc::dist(const Loc& loc1, const Loc& loc2) {
	if(isOverlap(loc1, loc2))
		return 0;
	else
		return loc1.start < loc2.start ? loc2.start - loc1.end + 1: loc1.start - loc2.end + 1;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
