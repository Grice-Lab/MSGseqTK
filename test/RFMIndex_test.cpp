/*
 * FMIndex_test.cpp
 *
 *  Created on: May 2, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include "RFMIndex.h"

using std::cout;
using std::endl;
using EGriceLab::MSGseqClean::DNAalphabet;
using EGriceLab::MSGseqClean::DNAseq;
using EGriceLab::MSGseqClean::RFMIndex;

int main() {
	cerr << "Max_length: " << RFMIndex::MAX_LENGTH << endl;
	DNAseq seqM1, seqM2;
	RFMIndex rfm1, rfm2;
	saidx_t count1 = 0, count2 = 0;
	DNAseq pat("CTAG");

	DNAseq seq1("CTAGCATAGAC");
	cout << "seq1:" << endl << seq1 << endl;
	seqM1 += seq1;
	rfm1 += RFMIndex(seq1);
	cout << "rfm1.length(): " << rfm1.length() << endl;
	cout << "rfm1.getBWT():" << endl << rfm1.getBWT() << endl;
	cout << "rfm1.getSeq():" << endl << rfm1.getSeq() << endl;
	if(rfm1.getSeq() != seqM1)
		return 1;
	count1 = rfm1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 1)
		return 1;

	DNAseq seq2("CTAGCATCGAC");
	cout << "seq2:" << endl << seq2 << endl;
	seqM1.push_back(0);
	seqM1 += seq2;
	rfm1 += RFMIndex(seq2);
	cout << "rfm1.length(): " << rfm1.length() << endl;
	cout << "rfm1.getBWT():" << endl << rfm1.getBWT() << endl;
	cout << "rfm1.getSeq():" << endl << rfm1.getSeq() << endl;
	if(rfm1.getSeq() != seqM1)
		return 1;
	count1 = rfm1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 2)
		return 1;

	DNAseq seq3("CTAGCATGGAC");
	cout << "seq3:" << endl << seq3 << endl;
	seqM1.push_back(0);
	seqM1 += seq3;
	rfm1 += RFMIndex(seq3);
	cout << "rfm1.length(): " << rfm1.length() << endl;
	cout << "rfm1.getBWT():" << endl << rfm1.getBWT() << endl;
	cout << "rfm1.getSeq():" << endl << rfm1.getSeq() << endl;
	if(rfm1.getSeq() != seqM1)
		return 1;
	count1 = rfm1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 3)
		return 1;

	DNAseq seq4("CTAGCATTGAC");
	cout << "seq4:" << endl << seq4 << endl;
	seqM1.push_back(0);
	seqM1 += seq4;
	rfm1 += RFMIndex(seq4);
	cout << "rfm1.length(): " << rfm1.length() << endl;
	cout << "rfm1.getBWT():" << endl << rfm1.getBWT() << endl;
	cout << "rfm1.getSeq():" << endl << rfm1.getSeq() << endl;
	if(rfm1.getSeq() != seqM1)
		return 1;
	count1 = rfm1.count(pat);
	cout << "found " << count1 << " of " << pat << " in " << seqM1 << endl;
	if(count1 != 4)
		return 1;

	DNAseq seq5("CTAGCAACTAG");
	cout << "seq5:" << endl << seq5 << endl;
	seqM2 += seq5;
	rfm2 += RFMIndex(seq5);
	cout << "rfm2.length(): " << rfm2.length() << endl;
	cout << "rfm2.getBWT():" << endl << rfm2.getBWT() << endl;
	cout << "rfm2.getSeq():" << endl << rfm2.getSeq() << endl;
	if(rfm2.getSeq() != seqM2)
		return 1;
	count2 = rfm2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 2)
		return 1;

	DNAseq seq6("CTAGCACCTAG");
	cout << "seq6:" << endl << seq6 << endl;
	seqM2.push_back(0);
	seqM2 += seq6;
	rfm2 += RFMIndex(seq6);
	cout << "rfm2.length(): " << rfm2.length() << endl;
	cout << "rfm2.getBWT():" << endl << rfm2.getBWT() << endl;
	cout << "rfm2.getSeq():" << endl << rfm2.getSeq() << endl;
	if(rfm2.getSeq() != seqM2)
		return 1;
	count2 = rfm2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 4)
		return 1;

	DNAseq seq7("CTAGCAGCTAG");
	cout << "seq7:" << endl << seq7 << endl;
	seqM2.push_back(0);
	seqM2 += seq7;
	rfm2 += RFMIndex(seq7);
	cout << "rfm2.length(): " << rfm2.length() << endl;
	cout << "rfm2.getBWT():" << endl << rfm2.getBWT() << endl;
	cout << "rfm2.getSeq():" << endl << rfm2.getSeq() << endl;
	if(rfm2.getSeq() != seqM2)
		return 1;
	count2 = rfm2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 6)
		return 1;

	DNAseq seq8("CTAGCATCTAG");
	cout << "seq8:" << endl << seq8 << endl;
	seqM2.push_back(0);
	seqM2 += seq8;
	rfm2 += RFMIndex(seq8);
	cout << "rfm2.length(): " << rfm2.length() << endl;
	cout << "rfm2.getBWT():" << endl << rfm2.getBWT() << endl;
	cout << "rfm2.getSeq():" << endl << rfm2.getSeq() << endl;
	if(rfm2.getSeq() != seqM2)
		return 1;
	count2 = rfm2.count(pat);
	cout << "found " << count2 << " of " << pat << " in " << seqM2 << endl;
	if(count2 != 8)
		return 1;

	DNAseq seqM = seqM1;
	seqM.push_back(0);
	seqM += seqM2;
	cout << "seqM:" << endl << seqM << endl;
	RFMIndex rfm = rfm1 + rfm2;
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	saidx_t count = rfm.count(pat);
	if(rfm.getSeq() != seqM)
		return 1;
	cout << "found " << count << " of " << pat << " in " << seqM << endl;
	if(count != count1 + count2)
		return 1;
}
