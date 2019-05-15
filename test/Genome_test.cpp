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
using namespace EGriceLab::MSGseqTK;

int main() {
	/* part 1, Genome test */
	DNAseq chr1_1 = dna::encode("ATCGNatcgnTCGANtcgan");
	DNAseq chr1_2 = dna::encode("TCGANtcganATCGNatcgn");
	Genome g1("g1", "Staphylococcus aureus");
	g1.addChrom("chr1_1", chr1_1);
	g1.addChrom("chr1_2", chr1_2);

	DNAseq chr2_1 = dna::encode("AAAAANgggggNCCCCCNtttttn");
	DNAseq chr2_2 = dna::encode("TTTTTNcccccNGGGGGNaaaaan");
	Genome g2("g2", "Homo sapiens");
	g2.addChrom("chr2_1", chr2_1);
	g2.addChrom("chr2_2", chr2_2);

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
		cerr << "g1.id: " << g1.getId() << " g1N.id: " << g1N.getId() << endl;
		cerr << "g1.name: " << g1.getName() << " g1N.name: " << g1N.getName() << endl;
		cerr << "g1.numChroms(): " << g1.numChroms() << endl;
		cerr << "g1N.numChroms(): " << g1N.numChroms() << endl;
		cerr << "g1.chr1 == g1N.chr1: " << (g1.getChrom(0) == g1N.getChrom(0)) << endl;
		cerr << "g1.chr1.name: " << g1.getChrom(0).name << endl;
		cerr << "g1N.chr1.name: " << g1N.getChrom(0).name << endl;
		cerr << "g1.chr1.name: " << g1.getChrom(0).name << endl;
		cerr << "g1.chr2 == g1N.chr2: " << (g1.getChrom(1) == g1N.getChrom(1)) << endl;
		return EXIT_FAILURE;
	}
	if(g2N != g2) {
		cerr << "loaded genome g2 differs from original copy" << endl;
		return EXIT_FAILURE;
	}

	cout << "Genome tests passed" << endl;

	/* part 2, MetaGenome test */
	// rename genomes
	out.str("");
	MetaGenome mtg1, mtg2;
	mtg1.addGenome(g1);
	mtg1.addGenome(g2);
	mtg2.addGenome(g1N);
	mtg2.addGenome(g2N);

	mtg1.updateIndex();
	cout << "mtg1 updated" << endl;

	mtg2.updateIndex();
	cout << "mtg2 updated" << endl;

	if(mtg1 != mtg2) {
		cerr << "Metagenome mtg1 and mtg2 don't match" << endl;
		return EXIT_FAILURE;
	}

	mtg2.getGenome(0).setId("g1N");
	mtg2.getGenome(1).setId("g2N");
	MetaGenome mtg = mtg1 + mtg2;
	if(mtg.BDSize() != mtg1.BDSize() + mtg1.BDSize()) {
		cerr << "Merged MetaGenome size doesn't match" << endl;
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

	cout << "MetaGenome tests passed" << endl;
}
