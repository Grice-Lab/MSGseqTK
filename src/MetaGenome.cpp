/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <sstream>
#include <boost/lexical_cast.hpp>
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
	return size;
}

size_t MetaGenome::numChroms() const {
	size_t N = 0;
	for(const Genome& genome : genomes)
		N += genome.numChroms();
	return N;
}

void MetaGenome::addGenome(const Genome& genome, const DNAseq& genomeSeq) {
	if(genome.size() != genomeSeq.length()) {
		warningLog << "genome size and seq length are different, ignore" << endl;
		return;
	}
	genomes.push_back(genome);
	seq += genomeSeq;
}

ostream& MetaGenome::save(ostream& out) const {
	/* save basic info */
	const size_t NG = numGenomes();
	out.write((const char*) &NG, sizeof(size_t));
	for(const Genome& genome : genomes)
		genome.save(out);
	/* save nt16 compressed seq */
	return seq.nt16Save(out);
}

istream& MetaGenome::load(istream& in, bool ignoreSeq) {
	/* load basic info */
	size_t NG = 0;
	in.read((char*) &NG, sizeof(size_t));
	genomes.resize(NG);
	for(size_t i = 0; i < NG; ++i)
		genomes[i].load(in);
	/* load seq in nt16 compressed */
	if(!ignoreSeq)
			seq.nt16Load(in);
	/* update index */
	updateIndex();
	return in;
}

MetaGenome& MetaGenome::operator+=(const MetaGenome& other) {
	genomes.insert(genomes.end(), other.genomes.begin(), other.genomes.end());
	updateIndex();
	return *this;
}

MetaGenome operator+(const MetaGenome& lhs, const MetaGenome& rhs) {
	MetaGenome mtg;
	mtg.genomes.insert(mtg.genomes.end(), lhs.genomes.begin(), lhs.genomes.end());
	mtg.genomes.insert(mtg.genomes.end(), rhs.genomes.begin(), rhs.genomes.end());
	mtg.updateIndex();
	return mtg;
}

void MetaGenome::updateIndex() {
	/* clear old data */
	genomeIds.clear();
	chromNames.clear();
	genomeId2Idx.clear();
	chromName2Idx.clear();
	chromIdx2GenomeIdx.clear();
	chromIdx2Nbefore.clear();
	genomeIdx2Loc.clear();
	chromIdx2Loc.clear();
	const size_t NG = numGenomes();
	const size_t NC = numChroms();
	size_t gid = 0;
	size_t cid = 0;
	uint64_t gStart = 0;
	genomeIds.reserve(NG);
	chromNames.reserve(NC);
	genomeIdx2Loc.reserve(NG);
	chromIdx2Loc.reserve(NC);
	chromIdx2GenomeIdx.reserve(NC);
	chromIdx2Nbefore.reserve(NC);
	for(Genome& genome : genomes) {
		uint64_t cStart = 0;
		size_t nid = 0; // # of chrom before in this genome
		if(genomeId2Idx.count(genome.id) > 0) { // genome.id must be unique
			cerr << "Error: Redundant genome " << genome.displayId() << " found in database" << endl;
			abort();
		}
		genomeIds.push_back(genome.id);
		genomeId2Idx[genome.id] = gid;
		for(Genome::Chrom& chr : genome.chroms) {
			if(chromName2Idx.count(chr.name)) {
				warningLog << "Non-unique chrom name '" << chr.name << "' from genome '" << genome.displayId() << "' found, ";
				chr.name = getChromId(genome.id, chr.name);
				warningLog << "replacing with '" << chr.name << "'" << endl;
			}
			chromNames.push_back(chr.name);
			chromName2Idx[chr.name] = cid;
			chromIdx2GenomeIdx.push_back(gid);
			chromIdx2Nbefore.push_back(nid);
			chromIdx2Loc.push_back(Loc(gStart + cStart, gStart + cStart + chr.size + 1)); // include the null terminal
			cid++;
			nid++;
			cStart += chr.size + 1;
		}
		assert(cStart == genome.size());
		genomeIdx2Loc.push_back(Loc(gStart, gStart + cStart));
		gid++;
		gStart += cStart;
	}
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
