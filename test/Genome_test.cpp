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
	Genome g1("g1", "Staphylococcus aureus");
	g1.addChrom(Genome::Chrom("chr1", DNAseq("ATCGNatcgnTCGANtcgan")));
	g1.addChrom(Genome::Chrom("chr2", DNAseq("TCGANtcganATCGNatcgn")));

	Genome g2("g2", "Homo sapiens");
	g2.addChrom(Genome::Chrom("chr1", DNAseq("AAAAANgggggNCCCCCNtttttn")));
	g2.addChrom(Genome::Chrom("chr2", DNAseq("TTTTTNcccccNGGGGGNaaaaan")));

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
		cerr << "g1.chr1.seq: " << g1.getChrom(0).seq << endl;
		cerr << "g1N.chr1.seq:" << g1N.getChrom(0).seq << endl;
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
	out.str("");
	MetaGenome mtg1, mtg2;
	mtg1.push_back(g2);
	mtg1.push_back(g1);
	mtg2.push_back(g2N);
	mtg2.push_back(g1N);
	mtg1.updateIndex();
	mtg2.updateIndex();

	if(mtg1 != mtg2) {
		cerr << "Metagenome mtg1 and mtg2 don't match" << endl;
		return EXIT_FAILURE;
	}

	size_t i1, i2;
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

	/* avoid genome ID collistions */
	mtg2.getGenome(0).setId("g1N");
	mtg2.getGenome(1).setId("g2N");
	MetaGenome mtg = mtg1 + mtg2;
	if(mtg.size() != mtg1.size() + mtg1.size()) {
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
