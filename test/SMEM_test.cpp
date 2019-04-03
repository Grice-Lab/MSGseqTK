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
	mtg.updateIndex();

	DNAseq genomeSeq = genome.getSeq();
	genomeSeq.pop_back();
	FMDIndex fmdidx = FMDIndex(genomeSeq, true);
	assert(mtg.size() == fmdidx.length());

	cout << "SMEM search ..." << endl;
	PrimarySeq read("ACGTAGTA", "seq1");
	SeedList seeds = SMEMS::findSeeds(&read, &mtg, &fmdidx, inf);
	for(const SeedPair& seed : seeds)
		if(!isValidSeed(genomeSeq, read.getSeq(), seed))
			return EXIT_FAILURE;
}

bool isValidSeed(const DNAseq& target, const DNAseq& query, const SeedPair& seed) {
	DNAseq tSeg = target.substr(seed.getStart(), seed.length());
	DNAseq qSeg = query.substr(seed.getFrom(), seed.length());
	if(seed.getStrand() == GLoc::REV)
		dna::revcom(qSeg);
	if(tSeg != qSeg) {
		cerr << "Unmatched SMEM seq db: " << tSeg << " read: " << qSeg << endl;
		return false;
	}
	return true;
}
