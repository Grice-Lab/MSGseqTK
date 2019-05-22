/*
 * FMDIndex_test.cpp
 *
 *  Created on: May 2, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include "FMDIndex.h"
#include "MetaGenome.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<GLoc>& locs);

int main() {
	/* test part 1, basic and merging function */
	cerr << "Max_length: " << FMDIndex::MAX_LENGTH << endl;
	DNAseq seqM1, seqM2;
	FMDIndex fmdidx1, fmdidx2;
	saidx_t count1 = 0, count2 = 0;
	DNAseq pat = dna::encode("GTAGC");
	vector<GLoc> locs;

	/* series operator+= tests */
	DNAseq seq1 = MetaGenome::getBDSeq(dna::encode("GTAGCATAGAG"));
	seqM1 += seq1;
	cout << "seqM1:" << endl << seqM1 << endl;
	fmdidx1 += FMDIndex(seq1, false);
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	assert(fmdidx1.getSeq() == seqM1);
	fmdidx1.buildSA();
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	locs = fmdidx1.locateAll(pat, GLoc::FWD);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq2 = MetaGenome::getBDSeq(dna::encode("GTAGCATCGAG"));
	seqM1 += seq2;
	cout << "seqM1:" << endl << seqM1 << endl;
	fmdidx1 += FMDIndex(seq2, false);
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	assert(fmdidx1.getSeq() == seqM1);
	fmdidx1.buildSA();
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	locs = fmdidx1.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq3 = MetaGenome::getBDSeq(dna::encode("GTAGCATGGAG"));
	seqM1 += seq3;
	cout << "seqM1:" << endl << seqM1 << endl;
	fmdidx1 += FMDIndex(seq3, false);
	fmdidx1.buildSA();
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	assert(fmdidx1.getSeq() == seqM1);
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;

	locs = fmdidx1.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq4 = MetaGenome::getBDSeq(dna::encode("GTAGCATTGAG"));
	seqM1 += seq4;
	cout << "seqM1:" << endl << seqM1 << endl;
	fmdidx1 += FMDIndex(seq4, false);
	fmdidx1.buildSA();
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	assert(fmdidx1.getSeq() == seqM1);
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	locs = fmdidx1.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	/* series operator+ tests */
	DNAseq seq5 = MetaGenome::getBDSeq(dna::encode("GTAGCAACTAG"));
	seqM2 = seq5 + seqM2;
	cout << "seqM2:" << endl << seqM2 << endl;
	fmdidx2 = FMDIndex(seq5, false) + fmdidx2;
	fmdidx2.buildSA();
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	assert(fmdidx2.getSeq() == seqM2);
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	locs = fmdidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq6 = MetaGenome::getBDSeq(dna::encode("GTAGCACCTAG"));
	seqM2 = seq6 + seqM2;
	cout << "seqM2:" << endl << seqM2 << endl;
	fmdidx2 = FMDIndex(seq6, false) + fmdidx2;
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	assert(fmdidx2.getSeq() == seqM2);
	fmdidx2.buildSA();
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	locs = fmdidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq7 = MetaGenome::getBDSeq(dna::encode("GTAGCAGCTAC"));
	seqM2 = seq7 + seqM2;
	cout << "seqM2:" << endl << seqM2 << endl;
	fmdidx2 = FMDIndex(seq7, false) + fmdidx2;
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	assert(fmdidx2.getSeq() == seqM2);
	fmdidx2.buildSA();
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	locs = fmdidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq8 = MetaGenome::getBDSeq(dna::encode("GTAGCATCTAC"));
	seqM2 = seq8 + seqM2;
	cout << "seqM2:" << endl << seqM2 << endl;
	fmdidx2 = FMDIndex(seq8, false) + fmdidx2;
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	assert(fmdidx2.getSeq() == seqM2);
	fmdidx2.buildSA();
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	locs = fmdidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	/* hierarchical merge test */
	DNAseq seqM = seqM1 + seqM2;
	cout << "seqM:" << endl << seqM << endl;
	FMDIndex fmdidx = fmdidx1 + fmdidx2;
	cout << "fmdidx.length(): " << fmdidx.length() << endl;
	cout << "fmdidx.getBWT():" << endl << fmdidx.getBWT() << endl;
	cout << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl;
	saidx_t count = fmdidx.count(pat);
	assert(fmdidx.getSeq() == seqM);
	fmdidx.buildSA();
	cout << "found " << count << " of " << pat << " in " << seqM << endl;
	if(count != count1 + count2)
		return EXIT_FAILURE;
	locs = fmdidx.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM, pat, locs))
		return EXIT_FAILURE;
	locs = fmdidx.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM, pat, locs))
		return EXIT_FAILURE;
	/* test part 2, IO function */
	ostringstream out;
	fmdidx.save(out);
	if(out.bad()) {
		cerr << "Failed to save fmdidx: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	istringstream in(out.str());
	FMDIndex fmdidxNew;
	fmdidxNew.load(in);
	if(in.bad()) {
		cerr << "Failed to load fmdidxNew: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(fmdidxNew.getBWT() != fmdidx.getBWT()) {
		cerr << "loaded fmdidx diffs from original copy" << endl;
		return EXIT_FAILURE;
	}
}

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<GLoc>& locs) {
	return std::all_of(locs.begin(), locs.end(),
			[&](const GLoc& loc) {
		return loc.getStrand() == GLoc::FWD && pat == seq.substr(loc.getStart(), loc.length()) ||
				loc.getStrand() == GLoc::REV && dna::revcom(pat) == seq.substr(loc.getStart(), loc.length());
	});
}
