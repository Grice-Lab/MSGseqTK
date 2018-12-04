/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <sstream>
#include "Loc.h"
#include "MetaGenome.h"
#include "StringUtils.h"
#include "ProgEnv.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::istringstream;
using std::pair;

uint64_t MetaGenome::size() const {
	uint64_t size = 0;
	for(const Genome& genome : genomes)
		size += genome.size();
	return size; /* include null terminal for each Genome */
}

size_t MetaGenome::numChroms() const {
	size_t N = 0;
	for(const Genome& genome : genomes)
		N += genome.numChroms();
	return N;
}

ostream& MetaGenome::save(ostream& out) const {
	/* save basic info */
	const size_t NG = numGenomes();
	out.write((const char*) &NG, sizeof(size_t));
	for(const deque<Genome>::value_type& genome : genomes)
		genome.save(out);

	return out;
}

istream& MetaGenome::load(istream& in) {
	/* load basic info */
	size_t NG = 0;
	in.read((char*) &NG, sizeof(size_t));
	genomes.resize(NG);
	for(size_t i = 0; i < NG; ++i)
		genomes[i].load(in);

	updateIndex();
	return in;
}

MetaGenome& MetaGenome::operator+=(const MetaGenome& other) {
	genomes.insert(genomes.end(), other.genomes.begin(), other.genomes.end());
	updateIndex();
	return *this;
}

void MetaGenome::updateIndex() {
	/* clear old data */
	genomeId2Idx.clear();
	chromName2Idx.clear();
	genomeIdx2Loc.clear();
	chromIdx2Loc.clear();
	size_t gid = 0;
	size_t cid = 0;
	uint64_t gStart = 0;
	for(Genome& genome : genomes) {
		uint64_t cStart = 0;
		if(genomeId2Idx.count(genome.id) != 0)
			cerr << "genome " << genome.id << " exists with idx: " << genomeId2Idx.at(genome.id) << endl;
		assert(genomeId2Idx.count(genome.id) == 0);
		genomeId2Idx[genome.id] = gid;
		for(Genome::Chrom& chr : genome.chroms) {
			if(chromName2Idx.count(chr.name)) {
				warningLog << "Redundant chrom name '" << chr.name << "' found in genome '" << genome.id << " replacing it with '" << genome.id + "." + chr.name << "'" << endl;
				chr.name = genome.id + "." + chr.name;
			}
			chromName2Idx[chr.name] = cid;
			chromIdx2Loc[cid] = Loc(cStart, cStart + chr.size + 1); // include the null terminal
			cid++;
			cStart += chr.size + 1;
		}
		assert(cStart == genome.size());
		genomeIdx2Loc[gid] = Loc(gStart, gStart + cStart);
		gid++;
		gStart += cStart;
	}
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
