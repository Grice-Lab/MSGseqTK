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
	DNAseq bdSeq;
	Genome g1("g1");
	g1.addChrom("chr1_1", chr1_1);
	g1.addChrom("chr1_2", chr1_2);
	MetaGenome mtg1;
	mtg1.addGenome(g1);
	mtg1.updateIndex();
	FMDIndex fmdidx1(mtg1.getBDSeq(), true);
	fmdidx1.clearBWT();
	bdSeq += mtg1.getBDSeq();
	cerr << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cerr << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl << "mtg1.getBDSeq():" << endl << mtg1.getBDSeq() << endl;
	assert(fmdidx1.getSeq() == bdSeq);

	const DNAseq chr2_1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr2_2 = chr2_1;
	Genome g2("g2");
	g2.addChrom("chr2_1", chr2_1);
	g2.addChrom("chr2_2", chr2_2);
	MetaGenome mtg2;
	mtg2.addGenome(g2);
	mtg2.updateIndex();
	FMDIndex fmdidx2(mtg2.getBDSeq(), true);
	fmdidx2.clearBWT();
	bdSeq += mtg2.getBDSeq();
	cerr << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cerr << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl << "mtg2.getBDSeq():" << endl << mtg2.getBDSeq() << endl;
	assert(fmdidx2.getSeq() == mtg2.getBDSeq());

	FMDIndex fmdidx = fmdidx1 + fmdidx2;
//	FMDIndex fmdidx = FMDIndex(mtg.getBDSeq(), true);
	MetaGenome mtg = mtg1 + mtg2;
	fmdidx.buildSA();
	fmdidx.clearBWT();
	cerr << "fmdidx.getBWT():" << endl << fmdidx.getBWT() << endl;
	cerr << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl << "mtg.getBDSeq():" << endl << mtg1.getBDSeq() + mtg2.getBDSeq() << endl;
	assert(fmdidx.getSeq() == mtg1.getBDSeq() + mtg2.getBDSeq());
	cerr << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl;
//	assert(mtg.BDSize() == fmdidx.length());

	/* simple SMEM/seed test */
	PrimarySeq read("ACGTAGTA", "seq1");
	SMEM_LIST smems = SMEM_LIST::findSMEMS(&read, &mtg, &fmdidx, inf);
	cout << "found " << smems.size() << " smems" << endl;
	for(const SMEM& smem : smems)
		for(const SeedPair& seed : smem.getSeeds())
			if(!isValidSeed(mtg.getSeq(seed.getTid()), read.getSeq(), seed))
				return EXIT_FAILURE;

	SeedList seeds = SMEM_LIST::findSeeds(&read, &mtg, &fmdidx, inf);
	cout << "found " << seeds.size() << " all smems seeds" << endl;
	for(const SeedPair& seed : seeds)
		if(!isValidSeed(mtg.getSeq(seed.getTid()), read.getSeq(), seed))
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
