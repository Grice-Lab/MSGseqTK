/*
 * MEM_test.cpp
 *
 *  Created on: Jun 13, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include "MEM.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

static bool isValidMEM(const DNAseq& db, const MEM& mem);

int main() {
	const DNAseq genomeDB("ACGTCGTAGTACTACGNACGTCGTAGTACTACG");

	FMIndex fmidx = FMIndex(genomeDB.reverse(), true);

	DNAseq read("ACGTAGTA");
	PrimarySeq seq(read, "seq1");

	MEM mem;
	int i = 0;
	while( mem.to < read.length() ) {
		mem = MEM::findMEM(&seq, &fmidx, i == 0 || !mem.empty() ? mem.to : mem.to + 1);
		mem.evaluate();
		mem.findLocs();
		cout << "mem " << i << " between db: " << genomeDB << " and read: " << read << " found at from: " << mem.from << " to: " << mem.to << endl;
		cout << "all matched locs:" << endl;
		for(const Loc& loc : mem.locs)
			cout << " " << loc;
		cout << endl;
		cout << "loglik: " << mem.loglik() << " likelihood: " << mem.pvalue()
				<< " evalue: " << mem.evalue() << endl;
		if(!isValidMEM(genomeDB, mem))
			return EXIT_FAILURE;
		i++;
	}

	cout << endl << "resetting read" << endl;
	mem.reset();
	read = "TACGNACGT";
	i = 0;
	while( mem.to < read.length() ) {
		mem = MEM::findMEM(&seq, &fmidx, i == 0 || !mem.empty() ? mem.to : mem.to + 1);
		mem.findLocs();
		cout << "N containing mem " << i << " between db: " << genomeDB << " and read: " << read << " found at from: " << mem.from << " to: " << mem.to << endl;
		cout << "all matched locs:" << endl;
		for(const Loc& loc : mem.locs)
			cout << " " << loc;
		cout << endl;
		cout << "loglik: " << mem.loglik() << " likelihood: " << mem.pvalue()
				<< " evalue: " << mem.evalue() << endl;
		if(!isValidMEM(genomeDB, mem))
			return EXIT_FAILURE;
		i++;
	}
}

bool isValidMEM(const DNAseq& db, const MEM& mem) {
	for(const Loc& loc : mem.locs) {
		DNAseq dbSeg = db.substr(loc.start, loc.length());
		DNAseq readSeg = mem.seq->getSeq().substr(mem.from, mem.length());
		if(dbSeg != readSeg) {
			cerr << "Unmatched MEM seq db: " << dbSeg << " read: " << readSeg << endl;
			return false;
		}
	}
	return true;
}
