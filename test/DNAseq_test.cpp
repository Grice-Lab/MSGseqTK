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
using std::cerr;
using std::endl;
using namespace EGriceLab::MSGseqTK;
using namespace EGriceLab::SAMtools;

int main() {
	DNAseq src1("ATCGNTCGANNNNNNNNNNNatcgntcga");
	cout << "src1:" << endl << src1 << endl;
	DNAseq dest1 = src1;
	dest1.removeGaps();
	if(dest1 != DNAseq("ATCGTCGAatcgtcga")) {
		cerr << "src1 gaps are not properly removed" << endl;
		return EXIT_FAILURE;
	}
	else
		cout << "src1.removeGaps():" << endl << dest1 << endl;

	const DNAseq src2("ATCGURYSWKMBDHVN");
	DNAseq dest2 = src2.revcom();
	if(dest2 != DNAseq("NBDHVKMWSRYACGAT")) {
		cerr << "src2 is not properly reversed-complemented" << endl;
		return EXIT_FAILURE;
	}
	else
		cout << "src2.revcom():" << endl << dest2 << endl;

	/* nt16 encode/decode tests */
	BAM::seq_str src1Nt16 = src1.nt16Encode();
	dest1 = DNAseq::nt16Decode(src1.length(), src1Nt16);
	if(dest1 != src1) {
		cerr << "nt16 encode/decode generated copy differs from original version" << endl <<
				"src1: " << src1 << endl <<
				"dest1:" << dest1 << endl;
		return EXIT_FAILURE;
	}
}
