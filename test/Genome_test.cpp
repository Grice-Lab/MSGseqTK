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
#include "MetaGenome.h"

using namespace std;
using namespace EGriceLab::MSGseqClean;

int main() {
	/* part 1, Genome test */
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
	cout << "Genome IO tests passed" << endl;

	/* part 2, MetaGenome test */
	out.str("");
	MetaGenome mtg1, mtg2;
	mtg1.push_front(g1);
	mtg1.push_front(g2);
	mtg2.push_front(g1N);
	mtg2.push_front(g2N);
	if(mtg1 != mtg2) {
		cerr << "Metagenome mtg1 and mtg2 don't match" << endl;
		return EXIT_FAILURE;
	}

	MetaGenome mtg = mtg1 + mtg2;
	if(mtg.getSize() != mtg1.getSize() + mtg1.getSize()) {
		cerr << "Merged MetaGenome size doesn't match" << endl;
		return EXIT_FAILURE;
	}
	if(mtg.getBaseCount() != mtg1.getBaseCount() + mtg2.getBaseCount()) {
		cerr << "Merged MetaGenome base count doesn't match" << endl;
		return EXIT_FAILURE;
	}

	mtg.save(out);
	if(out.bad()) {
		cerr << "Failed to save MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	in.str(out.str());
	MetaGenome mtgN;
	mtgN.load(in);
	if(in.bad()) {
		cerr << "Failed to load MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(mtgN != mtg) {
		cerr << "loaded MetaGenome mtg differs from original copy" << endl;
		return EXIT_FAILURE;
	}
	cout << "MetaGenome IO tests passed" << endl;
}
