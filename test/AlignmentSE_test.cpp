/*
 * AlignmentSE_test.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include "AlignmentSE.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;
using namespace EGriceLab::SAMtools;

int main(int argc, char* argv[]) {
	string prog = argv[0];
	/* construct default ScoreScheme */
	ScoreScheme ss;

	/* perfect match test */
	DNAseq query1("ATCG");
	DNAseq target1("AAAATCGCCC");
	cout << "querys and targets constructed" << endl;

	cout << "query1:  " << query1 << endl;
	cout << "target1: " << target1 << endl;
	/* build an alignment */
	AlignmentSE aln1(query1, target1, "r1", 0, &ss);
	cout << "AlignmentSE constructed" << endl;
	double maxScore1 = aln1.calculateScores();
	cout << "aln1 score calculated, maxScore: " << maxScore1 << endl;
	aln1.backTrace();
	cout << "aln1 back-traced " << endl;
	cout << "alnFrom: " << aln1.alnFrom << " alnTo: " << aln1.alnTo << " alnQLen: " << aln1.getAlnQLen() << endl;
	cout << "alnStart: " << aln1.alnStart << " alnEnd: " << aln1.alnEnd << " alnTLen: " << aln1.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln1.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln1.getAlnTSeq() << endl;
	cout << "aln1 cigar-str: " << BAM::decodeCigar(aln1.getAlnCigar()) << endl;
	cout << "aln1 MD:Z: " << aln1.getAlnMDTag() << endl;
	if(query1.length() != BAM::cigar2QLen(aln1.getAlnCigar())) {
		cerr << "query1 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln1.getAlnQLen() != AlignmentSE::mdTag2alnQLen(aln1.getAlnMDTag())) {
		cerr << "query1 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* perfect match in middle test */
	DNAseq query2("TTTATCGAAA");
	DNAseq target2("AAAATCGCCC");
	cout << "querys and targets constructed" << endl;

	cout << "query2:  " << query2 << endl;
	cout << "target2: " << target2 << endl;
	/* build an alignment */
	AlignmentSE aln2(query2, target2, "r2", 0, &ss);
	cout << "AlignmentSE constructed" << endl;
	double maxScore2 = aln2.calculateScores();
	cout << "aln2 score calculated, maxScore: " << maxScore2 << endl;
	aln2.backTrace();
	cout << "aln2 back-traced " << endl;
	cout << "alnFrom: " << aln2.alnFrom << " alnTo: " << aln2.alnTo << " alnQLen: " << aln2.getAlnQLen() << endl;
	cout << "alnStart: " << aln2.alnStart << " alnEnd: " << aln2.alnEnd << " alnTLen: " << aln2.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln2.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln2.getAlnTSeq() << endl;
	cout << "aln2 cigar-str: " << BAM::decodeCigar(aln2.getAlnCigar()) << endl;
	cout << "aln2 MD:Z: " << aln2.getAlnMDTag() << endl;
	if(query2.length() != BAM::cigar2QLen(aln2.getAlnCigar())) {
		cerr << "query2 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln2.getAlnQLen() != AlignmentSE::mdTag2alnQLen(aln2.getAlnMDTag())) {
		cerr << "query2 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* mis-match test */
	DNAseq query3("AAAATCGATCGATCGATCGTTT");
	DNAseq target3("GGGATCGATCAATCGATCGCCC");
	cout << "querys and targets constructed" << endl;

	cout << "query3:  " << query3 << endl;
	cout << "target3: " << target3 << endl;
	/* build an alignment */
	AlignmentSE aln3(query3, target3, "r3", 0, &ss);
	cout << "AlignmentSE constructed" << endl;
	double maxScore3 = aln3.calculateScores();
	cout << "aln3 score calculated, maxScore: " << maxScore3 << endl;
	aln3.backTrace();
	cout << "aln3 back-traced " << endl;
	cout << "alnFrom: " << aln3.alnFrom << " alnTo: " << aln3.alnTo << " alnQLen: " << aln3.getAlnQLen() << endl;
	cout << "alnStart: " << aln3.alnStart << " alnEnd: " << aln3.alnEnd << " alnTLen: " << aln3.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln3.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln3.getAlnTSeq() << endl;
	cout << "aln3 cigar-str: " << BAM::decodeCigar(aln3.getAlnCigar()) << endl;
	cout << "aln3 MD:Z: " << aln3.getAlnMDTag() << endl;
	if(query3.length() != BAM::cigar2QLen(aln3.getAlnCigar())) {
		cerr << "query3 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln3.getAlnQLen() != AlignmentSE::mdTag2alnQLen(aln3.getAlnMDTag())) {
		cerr << "query3 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* insertion match test */
	DNAseq  query4("AAAATCGATCAATCGATCGTTT");
	DNAseq target4("AAAATCGATCATCGATCGTTT");
	cout << "querys and targets constructed" << endl;

	cout << "query4:  " << query4 << endl;
	cout << "target4: " << target4 << endl;
	/* build an alignment */
	AlignmentSE aln4(query4, target4, "r4", 0, &ss);
	cout << "AlignmentSE constructed" << endl;
	double maxScore4 = aln4.calculateScores();
	cout << "aln4 score calculated, maxScore: " << maxScore4 << endl;

	aln4.backTrace();
	cout << "aln4 back-traced " << endl;
	cout << "alnFrom: " << aln4.alnFrom << " alnTo: " << aln4.alnTo << " alnQLen: " << aln4.getAlnQLen() << endl;
	cout << "alnStart: " << aln4.alnStart << " alnEnd: " << aln4.alnEnd << " alnTLen: " << aln4.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln4.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln4.getAlnTSeq() << endl;
	cout << "aln4 cigar-str: " << BAM::decodeCigar(aln4.getAlnCigar()) << endl;
	cout << "aln4 MD:Z: " << aln4.getAlnMDTag() << endl;
	if(query4.length() != BAM::cigar2QLen(aln4.getAlnCigar())) {
		cerr << "query4 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln4.getAlnQLen() != AlignmentSE::mdTag2alnQLen(aln4.getAlnMDTag())) {
		cerr << "query4 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}

	/* deletion match test */
	DNAseq  query5("AAAATCGATCGTCGATCGTTT");
	DNAseq target5("AAAATCGATCGATCGATCGTTT");
	cout << "querys and targets constructed" << endl;

	cout << "query5:  " << query5 << endl;
	cout << "target5: " << target5 << endl;
	/* build an alignment */
	AlignmentSE aln5(query5, target5, "r5", 0, &ss);
	cout << "AlignmentSE constructed" << endl;
	double maxScore5 = aln5.calculateScores();
	cout << "aln5 score calculated, maxScore: " << maxScore5 << endl;

	aln5.backTrace();
	cout << "aln5 back-traced " << endl;
	cout << "alnFrom: " << aln5.alnFrom << " alnTo: " << aln5.alnTo << " alnQLen: " << aln5.getAlnQLen() << endl;
	cout << "alnStart: " << aln5.alnStart << " alnEnd: " << aln5.alnEnd << " alnTLen: " << aln5.getAlnTLen() << endl;
	cout << "alnQSeq: " << aln5.getAlnQSeq() << endl;
	cout << "alnTSeq: " << aln5.getAlnTSeq() << endl;
	cout << "aln5 cigar-str: " << BAM::decodeCigar(aln5.getAlnCigar()) << endl;
	cout << "aln5 MD:Z: " << aln5.getAlnMDTag() << endl;
	if(query5.length() != BAM::cigar2QLen(aln5.getAlnCigar())) {
		cerr << "query5 length doesn't match cigar-str" << endl;
		return EXIT_FAILURE;
	}
	if(aln5.getAlnQLen() != AlignmentSE::mdTag2alnQLen(aln5.getAlnMDTag())) {
		cerr << "query5 aligned length doesn't match MD:Z tag" << endl;
		return EXIT_FAILURE;
	}
}
