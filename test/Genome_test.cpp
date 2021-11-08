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

	ostringstream gOut, sOut;
	g1.save(gOut);
	if(gOut.bad()) {
		cerr << "Failed to save genome '" << g1.getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	g1.saveSeq(sOut);
	if(sOut.bad()) {
		cerr << "Failed to save genome seq '" << g1.getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	g2.save(gOut);
	if(gOut.bad()) {
		cerr << "Failed to save genome '" << g2.getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	g2.saveSeq(sOut);
	if(sOut.bad()) {
		cerr << "Failed to save genome seq '" << g2.getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	Genome g1N, g2N;
	istringstream gIn(gOut.str());
	istringstream sIn(sOut.str());

	g1N.load(gIn);
	if(gIn.bad()) {
		cerr << "Failed to load genome g1N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	g1N.loadSeq(sIn);
	if(sIn.bad()) {
		cerr << "Failed to load genome seq g1N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	g2N.load(gIn);
	if(gIn.bad()) {
		cerr << "Failed to load genome g2N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	g2N.loadSeq(sIn);
	if(sIn.bad()) {
		cerr << "Failed to load genome seq g2N: " << ::strerror(errno) << endl;
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
	cout << "Genome tests passed" << endl;

	/* part 2, MetaGenome test */
	// rename genomes
	gOut.str("");
	sOut.str("");
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

	MetaGenome mtg = mtg1 + mtg2;
	cout << "mtg1 and mtg2 merged" << endl;
	/* check genome redundancy */
	size_t nGenome = mtg.renameRedundantGenomes();
	if(nGenome > 0)
		cout << "Renamed " << nGenome << " genome ids at metagenome level due to redundancy" << endl;

	size_t nChr = mtg.renameRedundantChroms();
	if(nChr > 0)
		cout << "Renamed " << nChr << " chromosome names at metagenome level due to redundancy" << endl;
	mtg.updateIndex();

	if(mtg.BDSize() != mtg1.BDSize() + mtg1.BDSize()) {
		cerr << "Merged MetaGenome size doesn't match" << endl;
		return EXIT_FAILURE;
	}

	mtg.save(gOut);
	if(gOut.bad()) {
		cerr << "Failed to save MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	mtg.saveSeq(sOut);
	if(sOut.bad()) {
		cerr << "Failed to save MetaGenome seq: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	gIn.str(gOut.str());
	sIn.str(sOut.str());
	MetaGenome mtgN;
	mtgN.load(gIn);
	if(gIn.bad()) {
		cerr << "Failed to load MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	mtgN.loadSeq(sIn);
	if(sIn.bad()) {
		cerr << "Failed to load MetaGenome seq: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(mtgN != mtg) {
		cerr << "loaded MetaGenome mtg differs from original copy" << endl;
		return EXIT_FAILURE;
	}

	cout << "MetaGenome tests passed" << endl;
}
