/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <sstream>
#include "MetaGenome.h"
#include "StringUtils.h"
#include "ProgEnv.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::istringstream;

uint64_t MetaGenome::size() const {
	uint64_t size = 0;
	for(const Genome& genome : genomes)
		size += genome.size();
	return size; /* include null terminal for each Genome */
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

uint64_t MetaGenome::getGenomeShift(size_t i) const {
	uint64_t shift = 0;
	for(size_t k = 0; k < i; ++k)
		shift += getGenome(k).size();
	return shift;
}

MetaGenome::GENOME_SHIFTMAP MetaGenome::getGenomeShift() const {
	GENOME_SHIFTMAP shiftMap;
	uint64_t shift = 0;
	for(size_t i = 0; i < numGenomes(); ++i) {
		shiftMap[i] = shift;
		shift += getGenome(i).size();
	}
	return shiftMap;
}

size_t MetaGenome::getLocId(uint64_t loc) const {
	uint64_t start = 0;
	size_t id = 0;
	for(const Genome genome: genomes) {
		for(const Genome::Chrom& chr : genome.getChroms()) {
			if(start <= loc && loc < start + chr.size)
				return id;
			start += chr.size + 1; // including null
		}
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
	genomes.resize(N);
	for(size_t i = 0; i < N; ++i)
		genomes[i].load(in);
	return in;
}

MetaGenome& MetaGenome::operator+=(const MetaGenome& other) {
	genomes.insert(genomes.end(), other.genomes.begin(), other.genomes.end());
	return *this;
}

size_t MetaGenome::countGenome(const string& genomeId) const {
	string gid = Genome::formatName(genomeId);
	size_t c = 0;
	for(const Genome& genome : genomes)
		if(genome.getId() == gid)
			c++;
	return c;
}

bool MetaGenome::hasGenome(const string& genomeId) const {
	string gid = Genome::formatName(genomeId);
	for(const Genome& genome : genomes)
		if(genome.getId() == gid)
			return true;
	return false;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
