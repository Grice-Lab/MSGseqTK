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
	FMDIndex fmdidx1(mtg1.getBDSeq());
	fmdidx1.buildSA();
	bdSeq += mtg1.getBDSeq();
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl << "mtg1.getBDSeq():" << endl << mtg1.getBDSeq() << endl;
	assert(fmdidx1.getSeq() == bdSeq);

	const DNAseq chr2_1 = dna::encode("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");
	const DNAseq chr2_2 = chr2_1;
	Genome g2("g2");
	g2.addChrom("chr2_1", chr2_1);
	g2.addChrom("chr2_2", chr2_2);
	MetaGenome mtg2;
	mtg2.addGenome(g2);
	mtg2.updateIndex();
	FMDIndex fmdidx2(mtg2.getBDSeq());
	fmdidx2.buildSA();
	bdSeq += mtg2.getBDSeq();
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl << "mtg2.getBDSeq():" << endl << mtg2.getBDSeq() << endl;
	assert(fmdidx2.getSeq() == mtg2.getBDSeq());

	FMDIndex fmdidx = fmdidx1 + fmdidx2;
//	FMDIndex fmdidx = FMDIndex(mtg.getBDSeq(), true);
	MetaGenome mtg = mtg1 + mtg2;
	mtg.updateIndex();
	fmdidx.buildSA();
	cout << "fmdidx.getBWT():" << endl << fmdidx.getBWT() << endl;
	cout << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl;
	cout << "mtg.getBDSeq():" << endl << mtg.getBDSeq() << endl;
	assert(fmdidx.getSeq() == mtg.getBDSeq());

	/* simple SMEM/seed test */
	PrimarySeq read("ACGTAGTA", "seq1");
	MEM_LIST mems = SMEM_LIST::findMEMS(&read, &mtg, &fmdidx, 1, inf);
	cout << "finding MEMS for " << read.getSeq() << endl;
	cout << "found " << mems.size() << " mems" << endl;
	for(const MEM& mem : mems) {
		cout << "mem: " << mem << endl;
		for(const SeedPair& seed : mem.getSeeds()) {
			if(!isValidSeed(mtg.getSeq(seed.getTid()), read.getSeq(), seed))
				return EXIT_FAILURE;
		}
	}

	/* all SMEM/seed test */
	SeedList seeds = SMEM_LIST::findSeeds(&read, &mtg, &fmdidx, 1, 0, inf);
	cout << "finding seeds for " << read.getSeq() << endl;
	cout << "found " << seeds.size() << " all smems seeds" << endl;
	for(const SeedPair& seed : seeds) {
		if(!isValidSeed(mtg.getSeq(seed.getTid()), read.getSeq(), seed))
			return EXIT_FAILURE;
	}
	cout << endl;
}

bool isValidSeed(const DNAseq& target, const DNAseq& query, const SeedPair& seed) {
	DNAseq tSeg = target.substr(seed.getStart(), seed.length());
	DNAseq qSeg = seed.getStrand() == GLoc::FWD ? query.substr(seed.getFrom(), seed.length()) :
			dna::revcom(query).substr(seed.getFrom(), seed.length());
	if(tSeg != qSeg) {
		cout << "Unmatched seed: " << seed << " between target:" << endl <<
				target << endl << "query:" << endl << query << endl <<
				"tSeg: " << tSeg << " qSeg: " << qSeg << endl;
		return false;
	}
	return true;
}
