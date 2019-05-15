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
	FMDIndex fmdidx;
	const DNAseq chr1_1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr1_2 = chr1_1;
	Genome g1("g1");
	g1.addChrom("chr1_1", chr1_1);
	g1.addChrom("chr1_2", chr1_2);
	MetaGenome mtg1;
	mtg1.addGenome(g1);
	mtg1.update();
	fmdidx += FMDIndex(mtg1.getBDSeq(), false);

	const DNAseq chr2_1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr2_2 = chr2_1;
	Genome g2("g2");
	g2.addChrom("chr2_1", chr2_1);
	g2.addChrom("chr2_2", chr2_2);
	MetaGenome mtg2;
	mtg2.addGenome(g2);
	mtg2.update();
	fmdidx += FMDIndex(mtg2.getBDSeq(), false);

	MetaGenome mtg = mtg1 + mtg2;
	cerr << "bdSeq:" << endl << mtg.getBDSeq() << endl;
//	FMDIndex fmdidx = FMDIndex(mtg.getBDSeq(), true);
	fmdidx.buildSA();

	const DNAseq& genomeSeq = mtg.getSeq();
	assert(mtg.BDSize() == fmdidx.length());

	PrimarySeq read("ACGTAGTA", "seq1");
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
		cerr << "Unmatched seed: " << seed << " between target:" << endl <<
				target << endl << "query:" << endl << query << endl <<
				"tSeg: " << tSeg << " qSeg: " << qSeg << endl;
		return false;
	}
	return true;
}
