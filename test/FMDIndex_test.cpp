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
	FMDIndex fmdidx1, fmdidx2;
	saidx_t count1 = 0, count2 = 0;
	DNAseq pat = dna::encode("CTAGC");
	vector<GLoc> locs;

	/* series operator+= tests */
	DNAseq seq1 = dna::encode("CTAGCATAGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATAGAC")) + DNAalphabet::GAP_BASE;
	cout << "seq1:" << endl << seq1 << endl;
	seqM1 += seq1;
	fmdidx1 += FMDIndex(seq1, false);
	fmdidx1.buildSA();
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmdidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	{
		FMDIndex fmd1(seqM1);
		if(fmdidx1.getBWT() != fmd1.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx1.getBWT() << endl << fmd1.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 1)
		return EXIT_FAILURE;
	{
		FMDIndex fmd1(seqM1);
		if(fmdidx1.getBWT() != fmd1.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx1.getBWT() << endl << fmd1.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
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

	DNAseq seq2 = dna::encode("CTAGCATCGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATCGAC")) + DNAalphabet::GAP_BASE;
	cout << "seq2:" << endl << seq2 << endl;
	seqM1 += seq2;
	fmdidx1 += FMDIndex(seq2, false);
	fmdidx1.buildSA();
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmdidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	{
		FMDIndex fmd1(seqM1);
		if(fmdidx1.getBWT() != fmd1.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx1.getBWT() << endl << fmd1.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 2)
		return EXIT_FAILURE;
	{
		FMDIndex fmd1(seqM1);
		if(fmdidx1.getBWT() != fmd1.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx1.getBWT() << endl << fmd1.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
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

	DNAseq seq3 = dna::encode("CTAGCATGGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATGGAC")) + DNAalphabet::GAP_BASE;
	cout << "seq3:" << endl << seq3 << endl;
	seqM1 += seq3;
	fmdidx1 += FMDIndex(seq3, false);
	fmdidx1.buildSA();
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmdidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	{
		FMDIndex fmd1(seqM1);
		if(fmdidx1.getBWT() != fmd1.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx1.getBWT() << endl << fmd1.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 3)
		return EXIT_FAILURE;

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

	DNAseq seq4 = dna::encode("CTAGCATTGAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATTGAC")) + DNAalphabet::GAP_BASE;
	cout << "seq4:" << endl << seq4 << endl;
	seqM1 += seq4;
	fmdidx1 += FMDIndex(seq4, false);
	fmdidx1.buildSA();
	cout << "fmdidx1.length(): " << fmdidx1.length() << endl;
	cout << "fmdidx1.getBWT():" << endl << fmdidx1.getBWT() << endl;
	cout << "fmdidx1.getSeq():" << endl << fmdidx1.getSeq() << endl;
	cout << "seqM1:" << endl << seqM1 << endl;
	if(fmdidx1.getSeq() != seqM1)
		return EXIT_FAILURE;
	{
		FMDIndex fmd1(seqM1);
		if(fmdidx1.getBWT() != fmd1.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx1.getBWT() << endl << fmd1.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count1 = fmdidx1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 4)
		return EXIT_FAILURE;
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
	DNAseq seq5 = dna::encode("CTAGCAACTAG") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCAACTAG")) + DNAalphabet::GAP_BASE;
	cout << "seq5:" << endl << seq5 << endl;
	seqM2 = seq5 + seqM2;
	fmdidx2 = FMDIndex(seq5, false) + fmdidx2;
	fmdidx2.buildSA();
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	{
		FMDIndex fmd2(seqM2);
		if(fmdidx2.getBWT() != fmd2.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx2.getBWT() << endl << fmd2.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 1)
		return EXIT_FAILURE;
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

	DNAseq seq6 = dna::encode("CTAGCACCTAG") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCACCTAG")) + DNAalphabet::GAP_BASE;
	cout << "seq6:" << endl << seq6 << endl;
	seqM2 = seq6 + seqM2;
	fmdidx2 = FMDIndex(seq6, false) + fmdidx2;
	fmdidx2.buildSA();
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	{
		FMDIndex fmd2(seqM2);
		if(fmdidx2.getBWT() != fmd2.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx2.getBWT() << endl << fmd2.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 2)
		return EXIT_FAILURE;
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

	DNAseq seq7 = dna::encode("CTAGCAGCTAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCAGCTAC")) + DNAalphabet::GAP_BASE;
	cout << "seq7:" << endl << seq7 << endl;
	seqM2 = seq7 + seqM2;
	fmdidx2 = FMDIndex(seq7, false) + fmdidx2;
	fmdidx2.buildSA();
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	{
		FMDIndex fmd2(seqM2);
		if(fmdidx2.getBWT() != fmd2.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx2.getBWT() << endl << fmd2.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 3)
		return EXIT_FAILURE;
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

	DNAseq seq8 = dna::encode("CTAGCATCTAC") + DNAalphabet::GAP_BASE + dna::revcom(dna::encode("CTAGCATCTAC")) + DNAalphabet::GAP_BASE;
	cout << "seq8:" << endl << seq8 << endl;
	seqM2 = seq8 + seqM2;
	fmdidx2 = FMDIndex(seq8, false) + fmdidx2;
	fmdidx2.buildSA();
	cout << "fmdidx2.length(): " << fmdidx2.length() << endl;
	cout << "fmdidx2.getBWT():" << endl << fmdidx2.getBWT() << endl;
	cout << "fmdidx2.getSeq():" << endl << fmdidx2.getSeq() << endl;
	cout << "seqM2:" << endl << seqM2 << endl;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	{
		FMDIndex fmd2(seqM2);
		if(fmdidx2.getBWT() != fmd2.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx2.getBWT() << endl << fmd2.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
	count2 = fmdidx2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 4)
		return EXIT_FAILURE;
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
	DNAseq seqM = seqM1 + seqM2 ;
	cout << "seqM:" << endl << seqM << endl;
	FMDIndex fmdidx = fmdidx1 + fmdidx2;
	fmdidx.buildSA();
	cout << "fmdidx.length(): " << fmdidx.length() << endl;
	cout << "fmdidx.getBWT():" << endl << fmdidx.getBWT() << endl;
	cout << "fmdidx.getSeq():" << endl << fmdidx.getSeq() << endl;
	saidx_t count = fmdidx.count(pat);
	if(fmdidx.getSeq() != seqM)
		return EXIT_FAILURE;
	if(fmdidx2.getSeq() != seqM2)
		return EXIT_FAILURE;
	{
		FMDIndex fmd(seqM);
		if(fmdidx.getBWT() != fmd.getBWT()) {
			std::cerr << "Unmatched merged BWT and single BWT" << std::endl <<
					fmdidx.getBWT() << endl << fmd.getBWT() << endl;
			return EXIT_FAILURE;
		}
	}
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
