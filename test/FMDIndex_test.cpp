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

bool isValidLocs(const DNAseq& seq, const DNAseq& pat, const vector<GLoc>& locs);

int main() {
	/* test part 1, basic and merging function */
	cerr << "Max_length: " << FMDIndex::MAX_LENGTH << endl;
	DNAseq seqM1, seqM2;
	FMDIndex fmidx1, fmidx2;
	saidx_t count1 = 0, count2 = 0;
	DNAseq pat("CTAGC");
	vector<GLoc> locs;

	/* series operator+= tests */
	DNAseq seq1 = DNAseq("CTAGCATAGAC") + DNAseq::DNAgap + DNAseq("CTAGCATAGAC").revcom();
	cout << "seq1:" << endl << seq1 << endl;
	seqM1 += seq1 + DNAseq::DNAgap;
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

	DNAseq seq2 = DNAseq("CTAGCATCGAC") + DNAseq::DNAgap + DNAseq("CTAGCATCGAC").revcom();
	cout << "seq2:" << endl << seq2 << endl;
	seqM1 += seq2 + DNAseq::DNAgap;
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

	DNAseq seq3 = DNAseq("CTAGCATGGAC") + DNAseq::DNAgap + DNAseq("CTAGCATGGAC").revcom();
	cout << "seq3:" << endl << seq3 << endl;
	seqM1 += seq3 + DNAseq::DNAgap;
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

	DNAseq seq4 = DNAseq("CTAGCATTGAC") + DNAseq::DNAgap + DNAseq("CTAGCATTGAC").revcom();
	cout << "seq4:" << endl << seq4 << endl;
	seqM1 += seq4 + DNAseq::DNAgap;
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
	DNAseq seq5 = DNAseq("CTAGCAACTAG") + DNAseq::DNAgap + DNAseq("CTAGCAACTAG").revcom();
	cout << "seq5:" << endl << seq5 << endl;
	seqM2 = seq5 + DNAseq::DNAgap + seqM2;
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

	DNAseq seq6 = DNAseq("CTAGCACCTAG") + DNAseq::DNAgap + DNAseq("CTAGCACCTAG").revcom();
	cout << "seq6:" << endl << seq6 << endl;
	seqM2 = seq6 + DNAseq::DNAgap + seqM2;
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

	DNAseq seq7 = DNAseq("CTAGCAGCTAC") + DNAseq::DNAgap + DNAseq("CTAGCAGCTAC").revcom();
	cout << "seq7:" << endl << seq7 << endl;
	seqM2 = seq7 + DNAseq::DNAgap + seqM2;
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

	DNAseq seq8 = DNAseq("CTAGCATCTAC") + DNAseq::DNAgap + DNAseq("CTAGCATCTAC").revcom();
	cout << "seq8:" << endl << seq8 << endl;
	seqM2 = seq8 + DNAseq::DNAgap + seqM2;
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
				loc.strand == GLoc::REV && pat.revcom() == seq.substr(loc.start, loc.length());
	});
}
