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

ostream& Genome::save(ostream& out) const {
	StringUtils::saveString(name, out);
	size_t NChrom = numChroms();
	out.write((const char*) &NChrom, sizeof(size_t));
	for(const chrmap_t::value_type& pair : chromSize) {
		StringUtils::saveString(pair.first, out);
		out.write((const char*) &pair.second, sizeof(chrmap_t::mapped_type));
	}
	return out;
}

istream& Genome::load(istream& in) {
	StringUtils::loadString(name, in);
	size_t NChrom = 0;
	in.read((char*) &NChrom, sizeof(size_t));
	for(size_t i = 0; i < NChrom; ++i) {
		string chr;
		StringUtils::loadString(chr, in);
		in.read((char*) &chromSize[chr], sizeof(chrmap_t::mapped_type)); /* chromSize[chr] will be default constructed first */
	}
	return in;
}

vector<string> Genome::getChroms() const {
	vector<string> chroms;
	chroms.reserve(chromSize.size());
	for(const chrmap_t::value_type& pair : chromSize)
		chroms.push_back(pair.first);
	return chroms;
}

uint64_t Genome::getSize() const {
	uint64_t size = 0;
	for(const chrmap_t::value_type& pair : chromSize)
		size += pair.second;
	return size;
}

bool operator==(const Genome& lhs, const Genome& rhs) {
	if(lhs.name != rhs.name)
		return false;
	if(lhs.chromSize.size() != rhs.chromSize.size())
		return false;

	for(const Genome::chrmap_t::value_type& pair : lhs.chromSize)
		if(!(rhs.chromSize.count(pair.first) != 0 && rhs.getChromSize(pair.first) == pair.second))
			return false;
	return true;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
