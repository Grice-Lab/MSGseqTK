/*
 * SeedChain_test.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "SMEM.h"
#include "SeedChain.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

static bool isValidChain(const DNAseq& target, const DNAseq& query, const SeedChain& chain);

int main() {
//	const DNAseq chr1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr1 = dna::encode("ACGTCGTAACGTGGTAACGTTGTA");
	const DNAseq chr2 = dna::encode("ACGTCGTAACGTGGTAACGTTGTA");
	Genome genome("db");
	genome.addChrom("chr1", chr1);
	genome.addChrom("chr2", chr2);
	MetaGenome mtg;
	mtg.addGenome(genome);
	mtg.updateIndex();

	const DNAseq& genomeSeq = mtg.getBDSeq();
	FMDIndex fmdidx = FMDIndex(mtg.getBDSeq()).buildSA();
	assert(mtg.BDSize() == fmdidx.length());

	PrimarySeq read("ACGTAGTA", "seq1");
	cout << "SMEM search between genome:" << endl << genomeSeq << endl << "seq1:" << endl << read.getSeq() << endl;
	SeedList seeds = SMEM_LIST::findSeeds(&read, &mtg, &fmdidx, 4, inf);
	cout << "found " << seeds.size() << " seeds:" << endl;
	ChainList chains = SeedChain::getChains(seeds, read.length() * 0.2, read.length() * 0.1);
	cout << "found " << chains.size() << " chains" << endl;
	for(const SeedChain& chain : chains) {
		cout << "chain size: " << chain.size() << " seeds: ";
		for(const SeedPair& seed : chain)
			cout << seed << " --> ";
		cout << endl;
		if(!isValidChain(genomeSeq, read.getSeq(), chain))
			return EXIT_FAILURE;
	}
}

bool isValidChain(const DNAseq& target, const DNAseq& query, const SeedChain& chain) {
	for(SeedChain::const_iterator seed = chain.begin(); seed != chain.end(); ++seed) {
		// check seed compatitability
		if(seed < chain.end() - 1 && !SeedPair::isCompatitable(*seed, *(seed+1))) {
			cerr << "Incompatitable seeds found between " << *seed << " and " << *(seed+1) << endl;
			return false;
		}
		// check seed matching
		DNAseq tSeg = target.substr(seed->getStart(), seed->length());
		DNAseq qSeg = seed->getStrand() == GLoc::FWD ? query.substr(seed->getFrom(), seed->length()) :
				dna::revcom(query).substr(seed->getFrom(), seed->length());
		if(tSeg != qSeg) {
			cerr << "Unmatched Seed between db: " << tSeg << " read: " << qSeg << endl;
			return false;
		}
	}
	return true;
}

