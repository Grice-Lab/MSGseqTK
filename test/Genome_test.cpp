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
	Genome g1("Staphylococcus aureus");
	g1.addChrom("chr1", DNAseq("ATCGNatcgnTCGANtcgan").length());
	g1.addChrom("chr2", DNAseq("TCGANtcganATCGNatcgn").length());

	Genome g2("Homo sapiens");
	g2.addChrom("chr1", DNAseq("AAAAANgggggNCCCCCNtttttn").length());
	g2.addChrom("chr2", DNAseq("TTTTTNcccccNGGGGGNaaaaan").length());

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

	size_t i1;
	if((i1 = g1.getChromIndex(0)) != 0) {
		cerr << "chromIndex of g1 at loc 0 should be 0 but found " << i1 << endl;
		return EXIT_FAILURE;
	}
	if((i1 = g1.getChromIndex(20)) != 0) {
		cerr << "chromIndex of g1 at loc 20 should be 0 but found " << i1 << endl;
		return EXIT_FAILURE;
	}
	if((i1 = g1.getChromIndex(21)) != 1) {
		cerr << "chromIndex of g1 at loc 21 should be 1 but found " << i1 << endl;
		return EXIT_FAILURE;
	}

	size_t i2;
	if((i2 = g2.getChromIndex(0)) != 0) {
		cerr << "chromIndex of g2 at loc 0 should be 0 but found " << i2 << endl;
		return EXIT_FAILURE;
	}
	if((i2 = g2.getChromIndex(24)) != 0) {
		cerr << "chromIndex of g2 at loc 24 should be 0 but found " << i2 << endl;
		return EXIT_FAILURE;
	}
	if((i2 = g2.getChromIndex(25)) != 1) {
		cerr << "chromIndex of g2 at loc 25 should be 1 but found " << i2 << endl;
		return EXIT_FAILURE;
	}

	cout << "Genome tests passed" << endl;

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

	size_t j1, j2;

	if((i1 = mtg1.getGenomeIndex(0)) != 0 && (j1 = mtg1.getChromIndex(0)) != 0) {
		cerr << "genomeIndex and chromIndex of mt1 at loc 0 should be 0, 0 but found " << i1 << ", " << j1 << endl;
		return EXIT_FAILURE;
	}
	if((i1 = mtg1.getGenomeIndex(25)) != 0 && (j1 = mtg1.getChromIndex(25)) != 1) {
		cerr << "genomeIndex and chromIndex of mt1 at loc 25 should be 0, 1 but found " << i1 << ", " << j1 << endl;
		return EXIT_FAILURE;
	}
	if((i1 = mtg1.getGenomeIndex(51)) != 1 && (j1 = mtg1.getChromIndex(51)) != 0) {
		cerr << "genomeIndex and chromIndex of mt1 at loc 51 should be 1, 0 but found " << i1 << ", " << j1 << endl;
		return EXIT_FAILURE;
	}
	if((i1 = mtg1.getGenomeIndex(71)) != 1 && (j1 = mtg1.getChromIndex(71)) != 0) {
		cerr << "genomeIndex and chromIndex of mt1 at loc 67 should be 1, 0 but found " << i1 << ", " << j1 << endl;
		return EXIT_FAILURE;
	}
	if((i1 = mtg1.getGenomeIndex(72)) != 1 && (j1 = mtg1.getChromIndex(72)) != 1) {
		cerr << "genomeIndex and chromIndex of mt1 at loc 72 should be 1, 1 but found " << i1 << ", " << j1 << endl;
		return EXIT_FAILURE;
	}

	MetaGenome mtg = mtg1 + mtg2;
	if(mtg.getSize() != mtg1.getSize() + mtg1.getSize()) {
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
