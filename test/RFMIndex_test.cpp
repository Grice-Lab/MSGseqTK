/*
 * FMIndex_test.cpp
 *
 *  Created on: May 2, 2018
 *      Author: zhengqi
 */

#include "../src/RFMIndex.h"

#include <iostream>

using std::cout;
using std::endl;
using EGriceLab::MSGseqClean::DNAalphabet;
using EGriceLab::MSGseqClean::DNAseq;
using EGriceLab::MSGseqClean::RFMIndex;

int main() {
	cerr << "Max_length: " << RFMIndex::MAX_LENGTH << endl;
	DNAseq seq1("CTAGCATAGAC");
	cout << "seq1:" << endl << seq1 << endl;

	RFMIndex rfm1(seq1);
	cout << "rfm1.length(): " << rfm1.length() << endl;
	cout << "rfm1.getBWT():" << endl << rfm1.getBWT() << endl;
	cout << "rfm1.getSeq():" << endl << rfm1.getSeq() << endl;
	if(rfm1.getSeq() != seq1)
		return 1;

	DNAseq pat("CTAG");
	saidx_t count1 = rfm1.count(pat);
	cout << "found " << count1 << " of " << pat << " in seq1 " << seq1 << endl;
	if(count1 != 1)
		return 1;

	DNAseq seq2("CTAGCATCGAC");
	cout << "seq2:" << endl << seq2 << endl;
	RFMIndex rfm2(seq2);
	cout << "rfm2.length(): " << rfm2.length() << endl;
	cout << "rfm2.getBWT():" << endl << rfm2.getBWT() << endl;
	cout << "rfm2.getSeq():" << endl << rfm2.getSeq() << endl;
	if(rfm2.getSeq() != seq2)
		return 1;

	saidx_t count2 = rfm2.count(pat);
	cout << "found " << count2 << " of " << pat << " in seq2 " << seq2 << endl;
	if(count2 != 1)
		return 1;

	RFMIndex rfm = rfm1 + rfm2;
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	DNAseq seq = seq1;
	seq.push_back(0);
	seq += seq2;
	if(rfm.getSeq() != seq)
		return 1;
	saidx_t count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << (seq1 + seq2) << endl;
	if(count != count1 + count2)
		return 1;

	DNAseq seq3("CTAGCATTGAC");
	cout << "seq3:" << endl << seq3 << endl;
	RFMIndex rfm3(seq3);
	if(rfm3.getSeq() != seq3)
		return 1;

	saidx_t count3 = rfm3.count(pat);
	cout << "found " << count3 << " of " << pat << " in seq3 " << seq3 << endl;
	if(count3 != 1)
		return 1;

	seq.push_back(0);
	seq += seq3;
	cout << "seq:" << endl << seq << endl;
	RFMIndex rfmMerged = RFMIndex(seq);
	cout << "rfmMerged.getBWT():" << endl << rfmMerged.getBWT() << endl;
	cout << "rfmMerged.getSeq():" << endl << rfmMerged.getSeq() << endl;

	rfm += rfm3;
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	if(rfm.getSeq() != seq)
		return 1;
	count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << seq << endl;
	if(count != count1 + count2 + count3)
		return 1;
}
