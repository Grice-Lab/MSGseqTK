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

static bool isValidSMEM(const DNAseq& db, const SMEM& smem);

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
	for(int64_t from = 0, to = 1; from < read.length(); from = to) {
		cout << "finding smem at from: " << from << endl;
		SMEMS smems = SMEM::findSMEMS(&read, &mtg, &fmdidx, from, to);
		SMEM::findLocs(smems);
		cout << "found " << smems.size() << " SMEMS between read: " << read.getSeq() << " and db:" << endl << fmdidx.getSeq() << endl << "from: " << from << " to:" << to << endl;
		for(const SMEM& smem : smems) {
			if(!isValidSMEM(genomeSeq, smem))
				return EXIT_FAILURE;
		}
	}
}

bool isValidSMEM(const DNAseq& db, const SMEM& smem) {
	for(const GLoc& loc : smem.locs) {
		DNAseq dbSeg = db.substr(loc.start, loc.length());
		DNAseq readSeg = smem.seq->getSeq().substr(smem.from, smem.length());
		if(loc.strand == GLoc::REV)
			dna::revcom(readSeg);
		if(dbSeg != readSeg) {
			cerr << "Unmatched SMEM seq db: " << dbSeg << " read: " << readSeg << endl;
			return false;
		}
	}
	return true;
}
