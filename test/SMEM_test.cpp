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
	const DNAseq chr1_1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr1_2 = chr1_1;
	Genome g1("g1");
	g1.addChrom("chr1_1", chr1_1);
	g1.addChrom("chr1_2", chr1_2);
	MetaGenome mtg1;
	mtg1.addGenome(g1);
	mtg1.updateIndex();
	FMDIndex fmdidx1(mtg1.getBDSeq(), true);
	fmdidx1.clearBWT();
	cerr << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cerr << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl << "mtg1.getBDSeq():" << endl << mtg1.getBDSeq() << endl;
	assert(fmdidx1.getSeq() == mtg1.getBDSeq());

	const DNAseq chr2_1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr2_2 = chr2_1;
	Genome g2("g2");
	g2.addChrom("chr2_1", chr2_1);
	g2.addChrom("chr2_2", chr2_2);
	MetaGenome mtg2;
	mtg2.addGenome(g2);
	mtg2.updateIndex();
	FMDIndex fmdidx2(mtg1.getBDSeq(), true);
	fmdidx2.clearBWT();
	cerr << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cerr << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl << "mtg2.getBDSeq():" << endl << mtg2.getBDSeq() << endl;
	assert(fmdidx2.getSeq() == mtg2.getBDSeq());

	MetaGenome mtg = mtg1 + mtg2;
	FMDIndex fmdidx = fmdidx1 + fmdidx2;
//	FMDIndex fmdidx = FMDIndex(mtg.getBDSeq(), true);
	fmdidx.buildSA();
	fmdidx.clearBWT();
	cerr << "fmdidx.getBWT():" << endl << fmdidx.getBWT() << endl;
	cerr << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl << "mtg.getBDSeq():" << endl << mtg.getBDSeq() << endl;
	assert(fmdidx.getSeq() == mtg.getBDSeq());

	const DNAseq& genomeSeq = mtg.getSeq();
	cerr << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl;
//	assert(mtg.BDSize() == fmdidx.length());

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
