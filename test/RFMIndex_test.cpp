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
	DNAseq seq;
	RFMIndex rfm;
	saidx_t count = 0;
	DNAseq pat("CTAG");

	DNAseq seq1("CTAGCATTGAC");
	cout << "seq1:" << endl << seq1 << endl;
	seq += seq1;
	rfm += RFMIndex(seq1);
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	if(rfm.getSeq() != seq)
		return 1;
	count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << seq << endl;
	if(count != 1)
		return 1;

	DNAseq seq2("CTAGCATAGAC");
	cout << "seq2:" << endl << seq2 << endl;
	seq.push_back(0);
	seq += seq2;
	rfm += RFMIndex(seq2);
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	if(rfm.getSeq() != seq)
		return 1;
	count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << seq << endl;
	if(count != 2)
		return 1;

	DNAseq seq3("CTAGCATGGAC");
	cout << "seq3:" << endl << seq3 << endl;
	seq.push_back(0);
	seq += seq3;
	rfm += RFMIndex(seq3);
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	if(rfm.getSeq() != seq)
		return 1;
	count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << seq << endl;
	if(count != 3)
		return 1;

	DNAseq seq4("CTAGCATCGAC");
	cout << "seq4:" << endl << seq4 << endl;
	seq.push_back(0);
	seq += seq4;
	rfm += RFMIndex(seq4);
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	if(rfm.getSeq() != seq)
		return 1;
	count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << seq << endl;
	if(count != 4)
		return 1;

	DNAseq seq5("CTAGCATCTAG");
	cout << "seq5:" << endl << seq5 << endl;
	seq.push_back(0);
	seq += seq5;
	rfm += RFMIndex(seq5);
	cout << "rfm.length(): " << rfm.length() << endl;
	cout << "rfm.getBWT():" << endl << rfm.getBWT() << endl;
	cout << "rfm.getSeq():" << endl << rfm.getSeq() << endl;
	if(rfm.getSeq() != seq)
		return 1;
	count = rfm.count(pat);
	cout << "found " << count << " of " << pat << " in " << seq << endl;
	if(count != 6)
		return 1;
}
