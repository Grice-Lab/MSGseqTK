/*
 * Alignment_test.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include "../src/Alignment.h"

#include <string>
#include <iostream>

using namespace std;
using namespace EGriceLab::MSGseqTK;
using namespace EGriceLab::SAMtools;

int main(int argc, char* argv[]) {
	string prog = argv[0];
	/* construct default ScoreScheme */
	ScoreScheme ss;

	/* perfect match test */
	const PrimarySeq r1("ATCG", "r1");
	const PrimarySeq rc1 = r1.revcom();

	DNAseq target1 = dna::encode("AAAATCGCCC");

	cout << "reads and targets constructed" << endl;

	cout << "read1:  " << r1.getSeq() << endl;
	cout << "target1: " << target1 << endl;
	/* build an alignment */
	Alignment aln1(&r1, &rc1, &target1, 0, GLoc::FWD,
			0, r1.length(), 0, target1.length(), &ss);
	cout << "Alignment constructed" << endl;
	aln1.calculateScores();
	cout << "aln1 score calculated, score: " << aln1.getScore() << endl;
	aln1.backTrace();
	cout << "aln1 back-traced " << endl;
	cout << "alnFrom: " << aln1.getAlnFrom() << " alnTo: " << aln1.getAlnTo() << " alnQLen: " << aln1.getAlnQLen() << endl;
	cout << "alnStart: " << aln1.getAlnStart() << " alnEnd: " << aln1.getAlnEnd() << " alnTLen: " << aln1.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln1.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln1.getAlnTSeq() << endl;
	cout << "aln1 cigar-str: " << BAM::decodeCigar(aln1.getAlnCigar()) << endl;
	cout << "aln1 MD:Z: " << aln1.getAlnMDTag() << endl;
	if(r1.length() != BAM::cigar2QLen(aln1.getAlnCigar())) {
		cerr << "read1 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln1.getAlnQLen() != Alignment::mdTag2alnQLen(aln1.getAlnMDTag())) {
		cerr << "read1 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* perfect match in middle test */
	const PrimarySeq r2("TTTATCGAAA", "r2");
	const PrimarySeq rc2 = r2.revcom();
	DNAseq target2 = dna::encode("AAAATCGCCC");

	cout << "reads and targets constructed" << endl;

	cout << "read2:  " << r2.getSeq() << endl;
	cout << "target2: " << target2 << endl;
	/* build an alignment */
	Alignment aln2(&r2, &rc2, &target2, 0, GLoc::FWD,
			0, r2.length(), 0, target2.length(), &ss);
	cout << "Alignment constructed" << endl;
	aln2.calculateScores();
	cout << "aln2 score calculated, score: " << aln2.getScore() << endl;
	aln2.backTrace();
	cout << "aln2 back-traced " << endl;
	cout << "alnFrom: " << aln2.getAlnFrom() << " alnTo: " << aln2.getAlnTo() << " alnQLen: " << aln2.getAlnQLen() << endl;
	cout << "alnStart: " << aln2.getAlnStart() << " alnEnd: " << aln2.getAlnEnd() << " alnTLen: " << aln2.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln2.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln2.getAlnTSeq() << endl;
	cout << "aln2 cigar-str: " << BAM::decodeCigar(aln2.getAlnCigar()) << endl;
	cout << "aln2 MD:Z: " << aln2.getAlnMDTag() << endl;
	if(r2.length() != BAM::cigar2QLen(aln2.getAlnCigar())) {
		cerr << "read2 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln2.getAlnQLen() != Alignment::mdTag2alnQLen(aln2.getAlnMDTag())) {
		cerr << "read2 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* mis-match test */
	const PrimarySeq r3("AAAATCGATCGATCGATCGTTT", "r3");
	const PrimarySeq rc3 = r3.revcom();
	DNAseq target3 = dna::encode("GGGATCGATCAATCGATCGCCC");
	cout << "reads and targets constructed" << endl;

	cout << "read3:  " << r3.getSeq() << endl;
	cout << "target3: " << target3 << endl;
	/* build an alignment */
	Alignment aln3(&r3, &rc3, &target3, 0, GLoc::FWD,
			0, r3.length(), 0, target3.length(), &ss);
	cout << "Alignment constructed" << endl;
	cout << "Alignment constructed" << endl;
	aln3.calculateScores();
	cout << "aln3 score calculated, score: " << aln3.getScore() << endl;
	aln3.backTrace();
	cout << "aln3 back-traced " << endl;
	cout << "alnFrom: " << aln3.getAlnFrom() << " alnTo: " << aln3.getAlnTo() << " alnQLen: " << aln3.getAlnQLen() << endl;
	cout << "alnStart: " << aln3.getAlnStart() << " alnEnd: " << aln3.getAlnEnd() << " alnTLen: " << aln3.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln3.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln3.getAlnTSeq() << endl;
	cout << "aln3 cigar-str: " << BAM::decodeCigar(aln3.getAlnCigar()) << endl;
	cout << "aln3 MD:Z: " << aln3.getAlnMDTag() << endl;
	if(r3.length() != BAM::cigar2QLen(aln3.getAlnCigar())) {
		cerr << "read3 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln3.getAlnQLen() != Alignment::mdTag2alnQLen(aln3.getAlnMDTag())) {
		cerr << "read3 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* insertion match test */
	const PrimarySeq r4("AAAATCGATCAATCGATCGTTT", "r4");
	const PrimarySeq rc4 = r4.revcom();
	DNAseq target4 = dna::encode("AAAATCGATCATCGATCGTTT");
	cout << "reads and targets constructed" << endl;

	cout << "read4:  " << r4.getSeq() << endl;
	cout << "target4: " << target4 << endl;
	/* build an alignment */
	Alignment aln4(&r4, &rc4, &target4, 0, GLoc::FWD,
			0, r4.length(), 0, target4.length(), &ss);
	cout << "Alignment constructed" << endl;
	aln4.calculateScores();
	cout << "aln4 score calculated, score: " << aln4.getScore() << endl;
	aln4.backTrace();
	cout << "aln4 back-traced " << endl;
	cout << "alnFrom: " << aln4.getAlnFrom() << " alnTo: " << aln4.getAlnTo() << " alnQLen: " << aln4.getAlnQLen() << endl;
	cout << "alnStart: " << aln4.getAlnStart() << " alnEnd: " << aln4.getAlnEnd() << " alnTLen: " << aln4.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln4.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln4.getAlnTSeq() << endl;
	cout << "aln4 cigar-str: " << BAM::decodeCigar(aln4.getAlnCigar()) << endl;
	cout << "aln4 MD:Z: " << aln4.getAlnMDTag() << endl;
	if(r4.length() != BAM::cigar2QLen(aln4.getAlnCigar())) {
		cerr << "read4 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln4.getAlnQLen() != Alignment::mdTag2alnQLen(aln4.getAlnMDTag())) {
		cerr << "read4 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* deletion match test */
	const PrimarySeq r5("AAAATCGATCGTCGATCGTTT", "r5");
	const PrimarySeq rc5 = r5.revcom();
	DNAseq target5 = dna::encode("AAAATCGATCGATCGATCGTTT");
	cout << "reads and targets constructed" << endl;

	cout << "read5:  " << r5.getSeq() << endl;
	cout << "target5: " << target5 << endl;
	/* build an alignment */
	Alignment aln5(&r5, &rc5, &target5, 0, GLoc::FWD,
			0, r5.length(), 0, target5.length(), &ss);
	cout << "Alignment constructed" << endl;
	aln5.calculateScores();
	cout << "aln5 score calculated, score: " << aln5.getScore() << endl;
	aln5.backTrace();
	cout << "aln5 back-traced " << endl;
	cout << "alnFrom: " << aln5.getAlnFrom() << " alnTo: " << aln5.getAlnTo() << " alnQLen: " << aln5.getAlnQLen() << endl;
	cout << "alnStart: " << aln5.getAlnStart() << " alnEnd: " << aln5.getAlnEnd() << " alnTLen: " << aln5.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln5.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln5.getAlnTSeq() << endl;
	cout << "aln5 cigar-str: " << BAM::decodeCigar(aln5.getAlnCigar()) << endl;
	cout << "aln5 MD:Z: " << aln5.getAlnMDTag() << endl;
	if(r5.length() != BAM::cigar2QLen(aln5.getAlnCigar())) {
		cerr << "read5 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln5.getAlnQLen() != Alignment::mdTag2alnQLen(aln5.getAlnMDTag())) {
		cerr << "read5 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}
}
