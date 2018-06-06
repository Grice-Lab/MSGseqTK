/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include "MetaGenome.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqClean {
using Eigen::Map;

uint64_t MetaGenome::getSize() const {
	uint64_t size = 0;
	for(const deque<Genome>::value_type & genome : genomes)
		size += genome.getSize();
	return size;
}

Vector4l MetaGenome::getBaseCount() const {
	Vector4l count = Vector4l::Zero();
	for(const deque<Genome>::value_type& genome : genomes)
		count += genome.getBaseCount();
	return count;
}

vector<Genome> MetaGenome::getGenomes() const {
	vector<Genome> genomes;
	std::copy(this->genomes.begin(), this->genomes.end(), genomes.begin());
	return genomes;
}

vector<string> MetaGenome::getGenomeNames() const {
	vector<string> names;
	for(const deque<Genome>::value_type genome : genomes)
		names.push_back(genome.getName());
	return names;
}

size_t MetaGenome::getGenomeIndex(uint64_t loc) const {
	uint64_t start = 0;
	for(deque<Genome>::const_iterator genome = genomes.begin(); genome != genomes.end(); ++genome) {
		if(start <= loc && loc <= start + genome->getSize() + 1) /* an additional null terminal included */
			return genome - genomes.begin();
		start += genome->getSize() + 1;
	}
	return -1;
}

ostream& MetaGenome::save(ostream& out) const {
	size_t N = numGenomes();
	out.write((const char*) &N, sizeof(size_t));
	for(const deque<Genome>::value_type& genome : genomes)
		genome.save(out);
	return out;
}

istream& MetaGenome::load(istream& in) {
	size_t N = 0;
	in.read((char*) &N, sizeof(size_t));
	for(size_t i = 0; i < N; ++i) {
		Genome genome;
		genome.load(in);
		push_back(genome);
	}
	return in;
}

MetaGenome& MetaGenome::operator+=(const MetaGenome& other) {
	genomes.insert(genomes.end(), other.genomes.begin(), other.genomes.end());
	return *this;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
