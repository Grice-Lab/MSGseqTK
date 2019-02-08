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
	out.write((const char*) &strand, sizeof(STRAND));
	return out;
}

istream& Loc::load(istream& in) {
	in.read((char*) &start, sizeof(int64_t));
	in.read((char*) &end, sizeof(int64_t));
	in.read((char*) &strand, sizeof(STRAND));
	return in;
}

ostream& Loc::write(ostream& out) const {
	out << start << '-' << end << ":" << decodeStrand(strand); /* start-end */
	return out;
}

istream& Loc::read(istream& in) {
	in >> start;
	in.ignore(1, '-');
	in >> end;
	in.ignore(1, ':');
	char s;
	in >> s; /* automatical ignore additional space */
	strand = encodeStrand(s);
	return in;
}

char Loc::decodeStrand(STRAND strand) {
	switch(strand) {
	case FWD:
		return '+';
	case REV:
		return '-';
	case UNK:
		return '.';
	default:
		return '.';
	}
}

Loc::STRAND Loc::encodeStrand(char s) {
	switch(s) {
	case '+':
		return FWD;
	case '-':
		return REV;
	case '.':
		return UNK;
	default:
		return UNK;
	}
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
