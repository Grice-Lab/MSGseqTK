/*
 * GenomeLoc.cpp
 *
 *  Created on: Feb 13, 2019
 *      Author: zhengqi
 */

#include "GLoc.h"

namespace EGriceLab {
namespace MSGseqTK {

ostream& GLoc::save(ostream& out) const {
	Loc::save(out);
	out.write((const char*) &tid, sizeof(int64_t));
	out.write((const char*) &strand, sizeof(STRAND));
	return out;
}

istream& GLoc::load(istream& in) {
	Loc::load(in);
	in.read((char*) &tid, sizeof(int64_t));
	in.read((char*) &strand, sizeof(STRAND));
	return in;
}

ostream& GLoc::write(ostream& out) const {
	out << tid << ":";
	Loc::write(out);
	out << ":" << decodeStrand(strand);
	return out;
}

istream& GLoc::read(istream& in) {
	in >> tid;
	in.ignore(1, ':');
	Loc::read(in);
	in.ignore(1, ':');
	char s;
	in >> s;
	strand = encodeStrand(s);
	return in;
}

char GLoc::decodeStrand(STRAND strand) {
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

GLoc::STRAND GLoc::encodeStrand(char s) {
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
