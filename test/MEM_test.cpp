/*
 * MEM_test.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "MEM.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;
using dna::operator<<;

static bool isValidMEM(const DNAseq& db, const MEM& mem);

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

	cout << "fwd MEM search ..." << endl;
	PrimarySeq read("ACGTAGTA", "seq1");
	MEM mem;
	for(int from = 0; from < read.length(); from = mem.to + 1) {
		mem = MEM::findMEMfwd(&read, &mtg, &fmdidx, from);
		mem.evaluate();
		mem.findLocs();
		cout << "mem between db: " << genomeSeq << " and read: " << read.getSeq() << " found at from: " << mem.from << " to: " << mem.to << endl;
		cout << "all matched locs:" << endl;
		for(const GLoc& loc : mem.locs)
			cout << " " << loc;
		cout << endl;
		cout << "loglik: " << mem.loglik() << " likelihood: " << mem.pvalue()
				<< " evalue: " << mem.evalue() << endl;
		if(!isValidMEM(genomeSeq, mem))
			return EXIT_FAILURE;
	}

	cout << "rev MEM search ..." << endl;
	read = PrimarySeq("TACGNACGT", "seq2");
	for(int to = read.length(); to > 0; to = mem.from - 1) {
		mem = MEM::findMEMrev(&read, &mtg, &fmdidx, to);
		mem.evaluate();
		mem.findLocs();
		cout << "N containing mem between db: " << genomeSeq << " and read: " << read.getSeq() << " found at from: " << mem.from << " to: " << mem.to << endl;
		cout << "all matched locs:" << endl;
		for(const GLoc& loc : mem.locs)
			cout << " " << loc;
		cout << endl;
		cout << "loglik: " << mem.loglik() << " likelihood: " << mem.pvalue()
				<< " evalue: " << mem.evalue() << endl;
		if(!isValidMEM(genomeSeq, mem))
			return EXIT_FAILURE;
	}
}

bool isValidMEM(const DNAseq& db, const MEM& mem) {
	for(const GLoc& loc : mem.locs) {
		DNAseq dbSeg = db.substr(loc.start, loc.length());
		DNAseq readSeg = mem.seq->getSeq().substr(mem.from, mem.length());
		if(loc.strand == GLoc::REV)
			dna::revcom(readSeg);
		if(dbSeg != readSeg) {
			cerr << "Unmatched MEM seq db: " << dbSeg << " read: " << readSeg << endl;
			return false;
		}
	}
	return true;
}
