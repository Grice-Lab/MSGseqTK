/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
#include "MetaGenome.h"
#include "StringUtils.h"
#include "ProgEnv.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {

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

ostream& MetaGenome::writeGFF(ostream& out, GFF::Version ver, const string& src) const {
	/* write each genome */
	for(const Genome& genome : genomes)
		genome.writeGFF(out, ver, src);
	return out;
}

ostream& MetaGenome::writeGFF(ostream& out, const GENOME_ANNOMAP& genomeAnnos, GFF::Version ver, const string& src) const {
	if(genomeAnnos.empty())
		return writeGFF(out, ver, src);
	/* write each genome with external annotations */
	for(const Genome& genome : genomes) {
		if(genomeAnnos.count(genome.getId())) { /* annotation exists */
			genome.writeGFF(out, genomeAnnos.at(genome.getId()), ver, src);
		}
		else {
			genome.writeGFF(out, ver, src);
		}
	}
	return out;
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
