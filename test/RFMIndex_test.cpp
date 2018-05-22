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
using EGriceLab::MSGseqClean::DNAseq;
using EGriceLab::MSGseqClean::RFMIndex;

int main() {
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
}
