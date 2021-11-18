/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <sstream>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include "MetaGenome.h"
#include "StringUtils.h"
#include "ProgEnv.h"
#include "MSGseqTKConst.h"
#include "MSGseqTK_main.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::istringstream;
using std::pair;
using std::unordered_map;
using std::unordered_set;

size_t MetaGenome::size() const {
	size_t size = 0;
	for(const Genome& genome : genomes)
		size += genome.size() + genome.numChroms(); // include null terminal
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
	return out;
}

ostream& MetaGenome::saveSeq(ostream& out) const {
	for(const Genome& genome : genomes)
		genome.saveSeq(out);
	return out;
}

istream& MetaGenome::load(istream& in) {
	/* load basic info */
	size_t NG = 0;
	in.read((char*) &NG, sizeof(size_t));
	genomes.resize(NG);
	for(size_t i = 0; i < NG; ++i)
		genomes[i].load(in);
	/* update index */
	updateIndex();
	return in;
}

istream& MetaGenome::loadSeq(istream& in) {
	for(Genome& genome : genomes)
		genome.loadSeq(in);
	return in;
}

DNAseq MetaGenome::loadSeq(size_t tid, istream& in) const {
	const int64_t tLen = getChrom(tid).len;
	DNAseq tSeq;
	in.seekg(getChromStart(tid)); // seek binary input to tid start
	StringUtils::loadString(tSeq, in, tLen);
	assert(0 == in.get());
	return tSeq;
}

DNAseq MetaGenome::loadBDSeq(size_t tStart, size_t tEnd, istream& in) const {
	assert(tStart < tEnd);
	DNAseq bdSeq(getChromBDLength(tStart, tEnd), DNAalphabet::GAP_BASE);
	DNAseq::iterator bdIt = bdSeq.begin();
	for(size_t tid = tStart; tid < tEnd; ++tid) {
		DNAseq seq = loadSeq(tid, in);
		const DNAseq::size_type N = seq.length();
		/* copy seq twice */
		std::copy(seq.begin(), seq.end(), bdIt);
		std::copy(seq.begin(), seq.end(), bdIt + N + 1);
		/* to basic */
		dna::toBasic(bdIt, bdIt + 2 * N + 1);
		/* revcom */
		dna::revcom(bdIt + N + 1, bdIt + 2 * N + 1);
		/* update */
		bdIt += 2 * N + 2;
	}
	return bdSeq;
}

DNAseq MetaGenome::getBDSeq(const DNAseq& seq) {
	const DNAseq::size_type N = seq.length();
	DNAseq bdSeq(2 * (N + 1), DNAalphabet::GAP_BASE); /* including gaps */
	/* copy seq twice */
	std::copy(seq.begin(), seq.end(), bdSeq.begin());
	std::copy(seq.begin(), seq.end(), bdSeq.begin() + N + 1);
	dna::toBasic(bdSeq); // use only non-ambiguous bases
	dna::revcom(bdSeq.begin() + N + 1, bdSeq.begin() + 2 * N + 1);
	return bdSeq;
}

MetaGenome& MetaGenome::append(const MetaGenome& other) {
	genomes.insert(genomes.end(), other.genomes.begin(), other.genomes.end());
	return *this;
}

MetaGenome& MetaGenome::prepend(const MetaGenome& other) {
	genomes.insert(genomes.begin(), other.genomes.begin(), other.genomes.end());
	return *this;
}

MetaGenome operator+(const MetaGenome& lhs, const MetaGenome& rhs) {
	MetaGenome mtgM(lhs);
	mtgM += rhs;
	return mtgM;
}

size_t MetaGenome::renameRedundantGenomes() {
	size_t n = 0;
	map<string, size_t> id2ver; /* id to version  map */
	for(Genome& genome : genomes) {
		if(id2ver.count(genome.id) > 0) { /* this genome id already exists */
			genome.id += "." + boost::lexical_cast<string>(id2ver[genome.id]);
			n++;
		}
		id2ver[genome.id]++;
	}
	return n;
}

size_t MetaGenome::removeRedundantGenomes() {
	size_t n = numGenomes();  /* old number of genomes */
	unordered_set<string> ids; /* unique genome id set */
	/* remove-erase redundnt genomes */
	genomes.erase(std::remove_if(genomes.begin(), genomes.end(),
			[&] (const Genome& genome)->bool { return !ids.insert(genome.id).second; }),
			genomes.end());
	return n - numGenomes(); /* return the difference */
}

size_t MetaGenome::renameRedundantChroms() {
	size_t n = 0;
	unordered_map<string, size_t> name2ver; /* id to version  map */
	for(Genome& genome : genomes) {
		for(Genome::Chrom& chr : genome.chroms) {
			if(name2ver.count(chr.name) > 0) { /* this genome id already exists */
				std::cerr << "found redundant chrom: " << chr.name << std::endl;
				chr.name += "." + boost::lexical_cast<string>(name2ver.at(chr.name));
				n++;
			}
			name2ver[chr.name]++;
		}
	}
	return n;
}

void MetaGenome::updateIndex() {
	/* clear old data */
	genomeIds.clear();
	chromNames.clear();
	genomeId2Idx.clear();
	chromName2Idx.clear();
	chromIdx2GenomeIdx.clear();
	chromIdx2Nbefore.clear();
//	genomeIdx2Loc.clear();
	chromIdx2Loc.clear();
	chromIdx2BDLoc.clear();
	const size_t NG = numGenomes();
	const size_t NC = numChroms();
	size_t gid = 0;
	size_t cid = 0;
	int64_t gStart = 0;
	int64_t gStart2 = 0;
	genomeIds.reserve(NG);
	chromNames.reserve(NC);
//	genomeIdx2Loc.reserve(NG);
	chromIdx2Loc.reserve(NC);
	chromIdx2BDLoc.reserve(NC);
	chromIdx2GenomeIdx.reserve(NC);
	chromIdx2Nbefore.reserve(NC);

	for(Genome& genome : genomes) {
		int64_t cStart = 0;
		int64_t cStart2 = 0;
		size_t nid = 0; // # of chrom before in this genome
		assert(genomeId2Idx.count(genome.id) == 0); // genome.id must be unique
		genomeIds.push_back(genome.id);
		genomeId2Idx[genome.id] = gid;
		for(Genome::Chrom& chr : genome.chroms) {
			assert(chromName2Idx.count(chr.name) == 0);
			chromNames.push_back(chr.name);
			chromName2Idx[chr.name] = cid;
			chromIdx2GenomeIdx.push_back(gid);
			chromIdx2Nbefore.push_back(nid);
			chromIdx2Loc.push_back(Loc(gStart + cStart, gStart + cStart + chr.size() + 1)); // include N terminal
			chromIdx2BDLoc.push_back(Loc(gStart2 + cStart2, gStart2 + cStart2 + (chr.size() + 1) * 2)); // include N terminal
			cid++;
			nid++;
			cStart += chr.size() + 1;
			cStart2 += (chr.size() + 1) * 2;
		}
		assert(cStart == genome.size() + genome.numChroms());
		assert(cStart2 == (genome.size() + genome.numChroms()) * 2);
//		genomeIdx2Loc.push_back(Loc(gStart, gStart + cStart));
		gid++;
		gStart += cStart;
		gStart2 += cStart2;
	}
	assert(gStart == size());
	assert(gStart2 == BDSize());
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
