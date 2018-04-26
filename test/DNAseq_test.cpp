/*
 * DNAseq-test.cpp
 *  Unit test for DNAseq class
 *  Created on: Apr 24, 2018
 *      Author: zhengqi
 */
#include <iostream>
#include <algorithm>
#include "DNAseq.h"

using std::cout;
using std::endl;
using EGriceLab::MSGSeqClean::DNAseq;

int main() {
	DNAseq src1("ATCGNTCGANNNNNNNNNNNatcgntcga");
	cout << "src1:" << endl << src1 << endl;
	DNAseq dest1 = src1;
	dest1.removeGaps();
	if(dest1 != DNAseq("ATCGTCGAatcgtcga"))
		return 1;
	else
		cout << "src1.removeGaps():" << endl << dest1 << endl;

	dest1 = src1;
	dest1.compressGaps();
	if(dest1.compressGaps() != DNAseq("ATCGNTCGANatcgntcga"))
		return 1;
	else
		cout << "src1.compressGaps():" << endl << dest1 << endl;

	DNAseq src2("ATCGURYSWKMBDHVN");
	if(src2 != DNAseq("ATCGTACGTGCTACGN"))
		return 1;
	else
		cout << "src2:" << endl << src2 << endl;

	DNAseq dest2 = src2.revcom();
	cout << "src2.revcom():" << endl << dest2 << endl;
	if(dest2 != DNAseq("NCGTAGCACGTACGAT"))
		return 1;
	else
		cout << "src2.revcom():" << endl << dest2 << endl;
}
