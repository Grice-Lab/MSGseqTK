/*
 * Loc.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: zhengqi
 */

#include "Loc.h"

namespace EGriceLab {
namespace MSGseqClean {

ostream& Loc::save(ostream& out) const {
	out.write((const char*) &start, sizeof(uint64_t));
	out.write((const char*) &end, sizeof(uint64_t));
	return out;
}

istream& Loc::load(istream& in) {
	in.read((char*) &start, sizeof(uint64_t));
	in.read((char*) &end, sizeof(uint64_t));
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

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
