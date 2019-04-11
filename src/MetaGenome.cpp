/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <sstream>
#include <utility>
#include <boost/lexical_cast.hpp>
#include "MetaGenome.h"
#include "StringUtils.h"
#include "ProgEnv.h"
#include "MSGseqTKConst.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::istringstream;
using std::pair;

size_t MetaGenome::size() const {
	size_t size = 0;
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

ostream& MetaGenome::save(ostream& out) const {
	/* save basic info */
	const size_t NG = numGenomes();
	out.write((const char*) &NG, sizeof(size_t));
	for(const Genome& genome : genomes)
		genome.save(out);

	/* save seqs */
	const size_t NS = seqs.size();
	out.write((const char*) &NS, sizeof(size_t));
	for(const DNAseq& seq : seqs)
		StringUtils::saveString(seq, out);
	return out;
}

istream& MetaGenome::load(istream& in) {
	/* load basic info */
	size_t NG = 0;
	in.read((char*) &NG, sizeof(size_t));
	genomes.resize(NG);
	for(size_t i = 0; i < NG; ++i)
		genomes[i].load(in);

	/* load seqs in reversed order */
	size_t NS = 0;
	in.read((char*) &NS, sizeof(size_t));
	assert(NS == numChroms());
	seqs.resize(NS);
	for(size_t i = 0; i < NS; ++i)
		StringUtils::loadString(seqs[i], in);

	/* update index */
	updateIndex();

	return in;
}

MetaGenome& MetaGenome::operator+=(const MetaGenome& other) {
	genomes.insert(genomes.end(), other.genomes.begin(), other.genomes.end());
	seqs.insert(seqs.end(), other.seqs.begin(), other.seqs.end());
	updateIndex();
	return *this;
}

MetaGenome operator+(const MetaGenome& lhs, const MetaGenome& rhs) {
	MetaGenome mtgM(lhs);
	mtgM += rhs;
	return mtgM;
}

void MetaGenome::updateSeq() {
	seqs.clear();
	for(Genome& genome : genomes) {
		for(Genome::Chrom& chr : genome.chroms) {
			seqs.push_back(std::move(chr.seq));
			chr.clearSeq();
		}
	}
}

void MetaGenome::updateIndex() {
	/* clear old data */
	genomeIds.clear();
	chromNames.clear();
	genomeId2Idx.clear();
	chromName2Idx.clear();
	chromIdx2GenomeIdx.clear();
	chromIdx2Nbefore.clear();
	chromIdx2endPos.clear();
//	genomeIdx2Loc.clear();
	chromIdx2Loc.clear();
	const size_t NS = seqs.size();
	const size_t NG = numGenomes();
	const size_t NC = numChroms();
	assert(NS == NC);
	size_t gid = 0;
	size_t cid = 0;
	int64_t gStart = 0;
	size_t endPos = 0; // distance to the end of the output/input
	size_t L = 0; // total length of sequences
	genomeIds.reserve(NG);
	chromNames.reserve(NC);
//	genomeIdx2Loc.reserve(NG);
	chromIdx2Loc.reserve(NC);
	chromIdx2GenomeIdx.reserve(NC);
	chromIdx2Nbefore.reserve(NC);
	chromIdx2endPos.resize(NS); // allocate now since it is accessed backward

	/* update genomeIdx2EndPos */
	for(size_t i = NS; i > 0; --i) {
		L += (seqs[i - 1].length() + 1) * 2;
		endPos += sizeof(size_t) + seqs[i - 1].length() * sizeof(DNAseq::value_type); // both length and sequences are saved by StringUtils
		chromIdx2endPos[i - 1] = endPos;
	}

	for(Genome& genome : genomes) {
		int64_t cStart = 0;
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
			chromIdx2Loc.push_back(Loc(gStart + cStart, gStart + cStart + (chr.size() + 1) * 2)); // use BDseq size
			cid++;
			nid++;
			cStart += (chr.size() + 1) * 2; // update cStart
		}
		assert(cStart == (genome.size() + genome.numChroms()) * 2);
//		genomeIdx2Loc.push_back(Loc(gStart, gStart + cStart));
		gid++;
		gStart += cStart;
	}
	assert(gStart == BDSize());
	assert(L == BDSize());
}

DNAseq MetaGenome::getBDSeq() const {
	DNAseq mtgSeq;
	mtgSeq.reserve(BDSize());
	for(const DNAseq& seq : seqs)
		mtgSeq += dna::toBasic(seq) + DNAalphabet::GAP_BASE + dna::revcom(dna::toBasic(seq)) + DNAalphabet::GAP_BASE;
	return mtgSeq;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
