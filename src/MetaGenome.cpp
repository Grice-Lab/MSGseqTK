/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "MetaGenome.h"
#include "StringUtils.h"
#include "ProgEnv.h"

namespace EGriceLab {
namespace MSGseqTK {
using Eigen::Map;

uint64_t MetaGenome::size() const {
	uint64_t size = 0;
	for(const Genome& genome : genomes)
		size += genome.size();
	return size; /* include null terminal for each Genome */
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
		if(start <= loc && loc < start + genome->size())
			return genome - genomes.begin();
		start += genome->size(); /* no null-terminal for genome */
	}
	return -1;
}

size_t MetaGenome::getChromIndex(uint64_t loc) const {
	size_t start = 0;
	for(const Genome& genome : genomes) {
		if(start <= loc && loc < start + genome.size())
			return genome.getChromIndex(loc - start);
		start += genome.size(); /* no null-terminal for genome */
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

bool operator==(const MetaGenome& lhs, const MetaGenome& rhs) {
	if(lhs.genomes.size() != rhs.genomes.size())
		return false;
	for(deque<Genome>::size_type i = 0; i < lhs.genomes.size(); ++i)
		if(lhs.genomes[i] != rhs.genomes[i])
			return false;
	return true;
}

ostream& MetaGenome::writeGFF(ostream& out, UCSC::GFF::Version ver, const string& src) const {
	/* write each genome */
	size_t shift = 0;
	for(const Genome& genome : genomes) {
		genome.writeGFF(out, ver, src, shift);
		shift += genome.size();
	}

	return out;
}

size_t MetaGenome::countGenome(const string& genomeName) const {
	string gname = Genome::formatName(genomeName);
	size_t c = 0;
	for(const Genome& genome : genomes)
		if(genome.getName() == gname)
			c++;
	return c;
}

bool MetaGenome::hasGenome(const string& genomeName) const {
	string gname = Genome::formatName(genomeName);
	for(const Genome& genome : genomes)
		if(genome.getName() == gname)
			return true;
	return false;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
