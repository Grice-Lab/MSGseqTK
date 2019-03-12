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

using namespace std;
using namespace EGriceLab::MSGseqTK;
using dna::operator<<;

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<GLoc>& locs);

int main() {
	/* test part 1, basic and merging function */
	cerr << "Max_length: " << FMDIndex::MAX_LENGTH << endl;
	DNAseq seqM1, seqM2;
	FMDIndex fmidx1, fmidx2;
	saidx_t count1 = 0, count2 = 0;
	DNAseq pat = dna::encode("CTAGC");
	vector<GLoc> locs;

	/* series operator+= tests */
	DNAseq seq1 = dna::encode("CTAGCATAGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATAGAC"));
	cout << "seq1:" << endl << seq1 << endl;
	seqM1 += seq1 + DNAalphabet::GAP_BASE;
	fmidx1 += FMDIndex(seq1, true);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 1)
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat, GLoc::FWD);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq2 = dna::encode("CTAGCATCGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATCGAC"));
	cout << "seq2:" << endl << seq2 << endl;
	seqM1 += seq2 + DNAalphabet::GAP_BASE;
	fmidx1 += FMDIndex(seq2);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 2)
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq3 = dna::encode("CTAGCATGGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATGGAC"));
	cout << "seq3:" << endl << seq3 << endl;
	seqM1 += seq3 + DNAalphabet::GAP_BASE;
	fmidx1 += FMDIndex(seq3);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 3)
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq4 = dna::encode("CTAGCATTGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATTGAC"));
	cout << "seq4:" << endl << seq4 << endl;
	seqM1 += seq4 + DNAalphabet::GAP_BASE;
	fmidx1 += FMDIndex(seq4);
	cout << "fmidx1.length(): " << fmidx1.length() << endl;
	cout << "fmidx1.getBWT():" << endl << fmidx1.getBWT() << endl;
	cout << "fmidx1.getSeq():" << endl << fmidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	count1 = fmidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 4)
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx1.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM1, pat, locs))
		return EXIT_FAILURE;

	/* series operator+ tests */
	DNAseq seq5 = dna::encode("CTAGCAACTAG") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCAACTAG"));
	cout << "seq5:" << endl << seq5 << endl;
	seqM2 = seq5 + DNAalphabet::GAP_BASE + seqM2;
	fmidx2 = FMDIndex(seq5, true) + fmidx2;
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 1)
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq6 = dna::encode("CTAGCACCTAG") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCACCTAG"));
	cout << "seq6:" << endl << seq6 << endl;
	seqM2 = seq6 + DNAalphabet::GAP_BASE + seqM2;
	fmidx2 = FMDIndex(seq6) + fmidx2;
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 2)
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq7 = dna::encode("CTAGCAGCTAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCAGCTAC"));
	cout << "seq7:" << endl << seq7 << endl;
	seqM2 = seq7 + DNAalphabet::GAP_BASE + seqM2;
	fmidx2 = FMDIndex(seq7) + fmidx2;
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 3)
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	DNAseq seq8 = dna::encode("CTAGCATCTAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATCTAC"));
	cout << "seq8:" << endl << seq8 << endl;
	seqM2 = seq8 + DNAalphabet::GAP_BASE + seqM2;
	fmidx2 = FMDIndex(seq8) + fmidx2;
	cout << "fmidx2.length(): " << fmidx2.length() << endl;
	cout << "fmidx2.getBWT():" << endl << fmidx2.getBWT() << endl;
	cout << "fmidx2.getSeq():" << endl << fmidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	count2 = fmidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 4)
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx2.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM2, pat, locs))
		return EXIT_FAILURE;

	/* hierarchical merge test */
	DNAseq seqM = seqM1 + seqM2 ;
	cout << "seqM:" << endl << seqM << endl;
	FMDIndex fmidx = fmidx1 + fmidx2;
	cout << "fmidx.length(): " << fmidx.length() << endl;
	cout << "fmidx.getBWT():" << endl << fmidx.getBWT() << endl;
	cout << "fmidx.getSeq():" << endl << fmidx.getSeq() << endl;
	saidx_t count = fmidx.count(pat);
	if(fmidx.getSeq() != seqM)
		return EXIT_FAILURE;
	cout << "found " << count << " of " << pat << " in " << seqM << endl;
	if(count != count1 + count2)
		return EXIT_FAILURE;
	locs = fmidx.locateAll(pat);
	cout << "All fwd loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM, pat, locs))
		return EXIT_FAILURE;
	locs = fmidx.locateAll(pat, GLoc::REV);
	cout << "All rev loc:";
	for(const GLoc& loc : locs)
		cout << " " << loc;
	cout << endl;
	if(!isValidLocs(seqM, pat, locs))
		return EXIT_FAILURE;
	/* test part 2, IO function */
	ostringstream out;
	fmidx.save(out);
	if(out.bad()) {
		cerr << "Failed to save fmidx: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	istringstream in(out.str());
	FMDIndex fmidxNew;
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

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<GLoc>& locs) {
	return std::all_of(locs.begin(), locs.end(),
			[&](const GLoc& loc) {
		return loc.strand == GLoc::FWD && pat == seq.substr(loc.start, loc.length()) ||
				loc.strand == GLoc::REV && dna::revcom(pat) == seq.substr(loc.start, loc.length());
	});
}
