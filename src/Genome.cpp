/*
 * Genome.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#include "Genome.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqClean {

ostream& Genome::Chrom::save(ostream& out) const {
	StringUtils::saveString(name, out);
	out.write((const char*) &size, sizeof(uint64_t));
	return out;
}

istream& Genome::Chrom::load(istream& in) {
	StringUtils::loadString(name, in);
	in.read((char*) &size, sizeof(uint64_t));
	return in;
}

bool Genome::hasChrom(const string& chrName) const {
	for(const Chrom& chr : chroms)
		if(chr.name == chrName)
			return true;
	return false;
}

uint64_t Genome::getChromSize(const string& chrName) const {
	for(const Chrom& chr : chroms)
		if(chr.name == chrName)
			return chr.size;
	return 0;
}

size_t Genome::getChromIndex(uint64_t loc) const {
	for(size_t i = 0, start = 0; i < chroms.size(); ++i, start += chroms[i].size + 1)
		if(i <= loc && loc <= i + chroms[i].size) /* include null terminal */
			return i;
	return -1;
}

ostream& Genome::save(ostream& out) const {
	StringUtils::saveString(name, out);
	size_t NChrom = numChroms();
	out.write((const char*) &NChrom, sizeof(size_t));
	for(const Chrom& chr : chroms)
		chr.save(out);
	return out;
}

istream& Genome::load(istream& in) {
	StringUtils::loadString(name, in);
	size_t NChrom = 0;
	in.read((char*) &NChrom, sizeof(size_t));
	for(size_t i = 0; i < NChrom; ++i) {
		chroms.push_back(Chrom()); /* default construct first */
		chroms[i].load(in);
	}
	return in;
}

uint64_t Genome::getSize() const {
	uint64_t size = 0;
	for(const Chrom& chr : chroms)
		size += chr.size;
	return size + numChroms(); /* add one null gap after each chromosome */
}

bool operator==(const Genome& lhs, const Genome& rhs) {
	if(lhs.name != rhs.name)
		return false;
	if(lhs.chroms.size() != rhs.chroms.size())
		return false;

	for(size_t i = 0; i < lhs.chroms.size(); ++i)
		if(lhs.chroms[i] != rhs.chroms[i])
			return false;
	return true;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
