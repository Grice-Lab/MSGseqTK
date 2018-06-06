/*
 * Genome_test.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include "Genome.h"

using namespace std;
using namespace EGriceLab::MSGseqClean;

int main() {
	Genome g1("Staphylococcus aureus", DNAseq("ATCGNatcgnTCGANtcga"));
	Genome g2("Homo sapiens", DNAseq("AAAAANgggggNCCCCCNttttt"));

	ostringstream out;
	g1.save(out);
	if(out.bad()) {
		cerr << "Failed to save '" << g1.getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	g2.save(out);
	if(out.bad()) {
		cerr << "Failed to save '" << g2.getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	Genome g1N, g2N;
	istringstream in(out.str());
	g1N.load(in);
	if(in.bad()) {
		cerr << "Failed to load g1N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	g2N.load(in);
	if(in.bad()) {
		cerr << "Failed to load g2N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(g1N != g1) {
		cerr << "loaded genome g1 differs from original copy" << endl;
		return EXIT_FAILURE;
	}
	if(g2N != g2) {
		cerr << "loaded genome g2 differs from original copy" << endl;
		return EXIT_FAILURE;
	}
}
