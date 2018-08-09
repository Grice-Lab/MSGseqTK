/*
 * GFF.cpp
 *
 *  Created on: Aug 7, 2018
 *      Author: zhengqi
 */

#include "StringUtils.h"
#include "GFF.h"

namespace EGriceLab {
namespace UCSC {

using namespace std;
const string GFF::GFF_SUFFIX = ".gff";
const string GFF::GTF_SUFFIX = ".gtf";
const string GFF::GFF3_SUFFIX = ".gff3";


istream& GFF::load(istream& in) {
	StringUtils::loadString(seqname, in);
	StringUtils::loadString(source, in);
	StringUtils::loadString(type, in);
	in.read((char *) &start, sizeof(long));
	in.read((char *) &end, sizeof(long));
	in.read((char *) &score, sizeof(double));
	in.read(&strand, 1);
	in.read((char *) &frame, sizeof(int));
	size_t nAttr = 0;
	string name, val;
	in.read((char *) &nAttr, sizeof(size_t));
	for(size_t i = 0; i < nAttr; ++i) {
		StringUtils::loadString(name, in);
		StringUtils::loadString(val, in);
		setAttr(name, val);
	}
	return in;
}

ostream& GFF::save(ostream& out) const {
	StringUtils::saveString(seqname, out);
	StringUtils::saveString(source, out);
	StringUtils::saveString(type, out);
	out.write((const char *) &start, sizeof(long));
	out.write((const char *) &end, sizeof(long));
	out.write((const char *) &score, sizeof(double));
	out.write(&strand, 1);
	out.write((const char *) &frame, sizeof(int));
	size_t nAttr = numAttrs();
	out.write((const char *) &nAttr, sizeof(size_t));
	for(attr_map::const_iterator pair = attrValues.begin(); pair != attrValues.end(); ++pair) {
		StringUtils::saveString(pair->first, out);
		StringUtils::saveString(pair->second, out);
	}
	return out;
}

istream& operator>>(istream& in, GFF& record) {
	std::getline(in, record.seqname, GFF::SEP);
	std::getline(in, record.source, GFF::SEP);
	std::getline(in, record.type, GFF::SEP);
	in >> record.start >> record.end >> record.score >> record.strand >> record.frame;
	string attrStr;
	std::getline(in, attrStr);
	record.readAttributes(attrStr);
	return in;
}

ostream& operator<<(ostream& out, const GFF& record) {
	out << record.seqname << GFF::SEP << record.source << GFF::SEP << record.type << GFF::SEP
		<< record.start << GFF::SEP << record.end << GFF::SEP << record.score << GFF::SEP
		<< record.strand << GFF::SEP << record.frame << GFF::SEP
		<< record.writeAttributes();
	return out;
}

GFF::VERSION GFF::guessVersion(const string& fn) {
	if(StringUtils::endsWith(fn, GTF_SUFFIX))
		return GTF;
	else if(StringUtils::endsWith(fn, GFF_SUFFIX) || StringUtils::endsWith(fn, GFF3_SUFFIX))
		return GFF3;
	else return UNK;
}

} /* namespace UCSC */
} /* namespace EGriceLab */

