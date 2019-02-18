/*
 * GenomeAnno_test.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include "GFF.h"
#include "Genome.h"
#include "MetaGenomeAnno.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

int main() {
	/* part 1, Genome test */
	Genome g1("Staphylococcus aureus");
	DNAseq chr1_1("ATCGNatcgnTCGANtcgan");
	DNAseq chr1_2("TCGANtcganATCGNatcgn");
	g1.addChrom("chr1", chr1_1);
	g1.addChrom("chr2", chr1_2);

	Genome g2("Homo sapiens");
	DNAseq chr2_1("AAAAANgggggNCCCCCNtttttn");
	DNAseq chr2_2("TTTTTNcccccNGGGGGNaaaaan");
	g2.addChrom("chr1", chr2_1);
	g2.addChrom("chr2", chr2_2);

	MetaGenomeAnno mta1;
	MetaGenomeAnno mta2;

	mta1.addGenome(g1);
	mta2.addGenome(g2);

	mta1.addChromAnno(MetaGenome::getChromId(g1.getId(), "chr1"), GFF(GFF::GFF3, "chr1", "test", "gene", 1, 8, 0, '.', 0));
	mta1.addChromAnno(MetaGenome::getChromId(g1.getId(), "chr1"), GFF(GFF::GFF3, "chr1", "test", "gene", 11, 18, 0, '.', 0));

	mta2.addChromAnno(MetaGenome::getChromId(g2.getId(), "chr1"), GFF(GFF::GFF3, "chr1", "test", "gene", 1, 8, 0, '.', 0));
	mta2.addChromAnno(MetaGenome::getChromId(g2.getId(), "chr1"), GFF(GFF::GFF3, "chr1", "test", "gene", 11, 18, 0, '.', 0));

	ostringstream out;
	mta1.save(out);
	if(out.bad()) {
		cerr << "Failed to save mta1: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	mta2.save(out);
	if(out.bad()) {
		cerr << "Failed to save mta2: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* binary IO tests */
	MetaGenomeAnno mta1N, mta2N;
	istringstream in(out.str());
	mta1N.load(in);
	if(in.bad()) {
		cerr << "Failed to load a1N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	mta2N.load(in);
	if(in.bad()) {
		cerr << "Failed to load a2N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(mta1 != mta1N) {
		cerr << "Loaded mta1 doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}
	if(mta2 != mta2N) {
		cerr << "Loaded mta1 doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}

	MetaGenomeAnno mta = mta1 + mta2;
	mta.save(out);
	if(out.bad()) {
		cerr << "Failed to save mta: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	in.str(out.str());
	MetaGenomeAnno mtaN;
	mtaN.load(in);
	if(in.bad()) {
		cerr << "Failed to load mtaN: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(mtaN != mta) {
		cerr << "Loaded mtaN doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}

	cout << "MetaGenomeAnno tests passed" << endl;
}
