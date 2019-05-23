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

int main() {
	const DNAseq src = dna::encode("ATCGURYSWKMBDHVN");
	DNAseq dest = dna::revcom(src);
	if(dest != dna::encode("NBDHVKMWSRYACGAT")) {
		cerr << "src is not properly reversed-complemented" << endl;
		return EXIT_FAILURE;
	}
	else
		cout << "src.revcom():" << endl << dest << endl;

	/* nt16 encode/decode tests */
	DNAseq srcNt16 = dna::nt16Encode(src);
	dest = dna::nt16Decode(src.length(), srcNt16);
	if(dest != src) {
		cerr << "nt16 encode/decode generated copy differs from original version" << endl <<
				"src: " << src << endl <<
				"dest:" << dest << endl;
		return EXIT_FAILURE;
	}
}
