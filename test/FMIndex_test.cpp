/*
 * FMIndex_test.cpp
 *
 *  Created on: May 2, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include "FMIndex.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<Loc>& locs);

int main() {
	/* test part 1, basic and merging function */
	cerr << "Max_length: " << FMIndex::MAX_LENGTH << endl;
	DNAseq seqM1, seqM2;
	FMIndex fmidx1, fmidx2;
	saidx_t count1 = 0, count2 = 0;
	DNAseq pat("CTAG");
	vector<Loc> locs;

	DNAseq seq1("CTAGCATAGAC");
	cout << "seq1:" << endl << seq1 << endl;
	seqM1 += seq1;
	fmidx1 += FMIndex(seq1, true);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 1)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx1.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, fmidx1.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq2("CTAGCATCGAC");
	cout << "seq2:" << endl << seq2 << endl;
	seqM1.push_back(0);
	seqM1 += seq2;
	fmidx1 += FMIndex(seq2);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 2)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx1.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, fmidx1.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq3("CTAGCATGGAC");
	cout << "seq3:" << endl << seq3 << endl;
	seqM1.push_back(0);
	seqM1 += seq3;
	fmidx1 += FMIndex(seq3);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 3)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx1.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, fmidx1.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq4("CTAGCATTGAC");
	cout << "seq4:" << endl << seq4 << endl;
	seqM1.push_back(0);
	seqM1 += seq4;
	fmidx1 += FMIndex(seq4);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 4)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx1.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, fmidx1.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq5("CTAGCAACTAG");
	cout << "seq5:" << endl << seq5 << endl;
	seqM2 += seq5;
	fmidx2 += FMIndex(seq5, true);
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 2)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx2.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, fmidx2.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq6("CTAGCACCTAG");
	cout << "seq6:" << endl << seq6 << endl;
	seqM2.push_back(0);
	seqM2 += seq6;
	fmidx2 += FMIndex(seq6);
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 4)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx2.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, fmidx2.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq7("CTAGCAGCTAG");
	cout << "seq7:" << endl << seq7 << endl;
	seqM2.push_back(0);
	seqM2 += seq7;
	fmidx2 += FMIndex(seq7);
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 6)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx2.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, fmidx2.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seq8("CTAGCATCTAG");
	cout << "seq8:" << endl << seq8 << endl;
	seqM2.push_back(0);
	seqM2 += seq8;
	fmidx2 += FMIndex(seq8);
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 8)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx2.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, fmidx2.locateAll(pat)))
		return EXIT_FAILURE;

	DNAseq seqM = seqM1;
	seqM.push_back(0);
	seqM += seqM2;
	cout << "seqM:" << endl << seqM << endl;
	FMIndex fmidx = fmidx1 + fmidx2;
	cout << "fmidx.length(): " << fmidx.length() << endl;
	cout << "fmidx.getBWT():" << endl << fmidx.getBWT() << endl;
	cout << "fmidx.getSeq():" << endl << fmidx.getSeq() << endl;
	saidx_t count = fmidx.count(pat);
	if(fmidx.getSeq() != seqM)
		return EXIT_FAILURE;
	cout << "found " << count << " of " << pat << " in " << seqM << endl;
	if(count != count1 + count2)
		return EXIT_FAILURE;
	cout << "All mapped loc:";
	for(const Loc& loc : fmidx.locateAll(pat))
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM, pat, fmidx.locateAll(pat)))
		return EXIT_FAILURE;

	/* test part 2, IO function */
	ostringstream out;
	fmidx.save(out);
	if(out.bad()) {
		cerr << "Failed to save fmidx: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	istringstream in(out.str());
	FMIndex fmidxNew;
	fmidxNew.load(in);
	if(in.bad()) {
		cerr << "Failed to load fmidxNew: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(fmidxNew.getBWTStr() != fmidx.getBWTStr()) {
		cerr << "loaded fmidx diffs from original copy" << endl;
		return EXIT_FAILURE;
	}
}

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<Loc>& locs) {
	for(const Loc& loc : locs)
		if(seq.substr(loc.start, loc.length()) != pat)
			return false;
	return true;
}
