/*
 * SSMEMS_test.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "SMEM.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

static bool isValidSeed(const DNAseq& target, const DNAseq& query, const SeedPair& seed);

int main() {
	const DNAseq chr1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr2 = chr1;
	Genome genome("db");
	genome.addChrom("chr1", chr1);
	genome.addChrom("chr2", chr2);
	MetaGenome mtg;
	mtg.addGenome(genome);
	mtg.update();

	const DNAseq& genomeSeq = mtg.getSeq();
	FMDIndex fmdidx = FMDIndex(mtg.getBDSeq(), true);
	cerr << "mtg.BDSize():" << mtg.BDSize() << " fmdidx.length():" << fmdidx.length() << endl;
	assert(mtg.BDSize() == fmdidx.length());

	PrimarySeq read("ACGTAGTA", "seq1");
	cout << "SMEM search between genome:" << endl << genomeSeq << endl <<
			"bdSeq:" << endl << mtg.getBDSeq() << endl << "query:" << endl << read.getSeq() << endl;
	SMEMS smems = SMEMS::findSMEMS(&read, &mtg, &fmdidx, inf);
	cout << "found " << smems.size() << " smems" << endl;
	for(const SMEM& smem : smems)
		for(const SeedPair& seed : smem.getSeeds())
			if(!isValidSeed(genomeSeq, read.getSeq(), seed))
				return EXIT_FAILURE;

	SeedList seeds = SMEMS::findSeeds(&read, &mtg, &fmdidx, inf);
	cout << "found " << seeds.size() << " all smems seeds" << endl;
	for(const SeedPair& seed : seeds)
		if(!isValidSeed(genomeSeq, read.getSeq(), seed))
			return EXIT_FAILURE;
	cout << endl;
}

bool isValidSeed(const DNAseq& target, const DNAseq& query, const SeedPair& seed) {
	DNAseq tSeg = target.substr(seed.getStart(), seed.length());
	DNAseq qSeg = seed.getStrand() == GLoc::FWD ? query.substr(seed.getFrom(), seed.length()) :
			dna::revcom(query).substr(seed.getFrom(), seed.length());
	if(tSeg != qSeg) {
		cerr << "Unmatched SMEM seq db: " << tSeg << " read: " << qSeg << endl;
		return false;
	}
	return true;
}
