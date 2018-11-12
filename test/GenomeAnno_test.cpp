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
	g1.addChrom("chr1", DNAseq("ATCGNatcgnTCGANtcgan").length());
	g1.addChrom("chr2", DNAseq("TCGANtcganATCGNatcgn").length());

	Genome g2("Homo sapiens");
	g2.addChrom("chr1", DNAseq("AAAAANgggggNCCCCCNtttttn").length());
	g2.addChrom("chr2", DNAseq("TTTTTNcccccNGGGGGNaaaaan").length());

	GenomeAnno a1(g1);
	a1.addAnno("chr1", GFF(GFF::GFF3, "chr1", "test", "gene", 1, 8, 0, '.', 0));
	a1.addAnno("chr1", GFF(GFF::GFF3, "chr1", "test", "gene", 11, 18, 0, '.', 0));

	GenomeAnno a2(g2);
	a2.addAnno("chr1", GFF(GFF::GFF3, "chr1", "test", "gene", 1, 8, 0, '.', 0));
	a2.addAnno("chr1", GFF(GFF::GFF3, "chr1", "test", "gene", 11, 18, 0, '.', 0));

	ostringstream out;
	a1.save(out);
	if(out.bad()) {
		cerr << "Failed to save '" << a1.getGenome().getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	a2.save(out);
	if(out.bad()) {
		cerr << "Failed to save '" << a2.getGenome().getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* binary IO tests */
	GenomeAnno a1N, a2N;
	istringstream in(out.str());
	a1N.load(in);
	if(in.bad()) {
		cerr << "Failed to load a1N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	a2N.load(in);
	if(in.bad()) {
		cerr << "Failed to load a2N: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	ostringstream outN;
	a1N.save(outN);
	if(outN.bad()) {
		cerr << "Failed to save '" << a1N.getGenome().getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	a2N.save(outN);
	if(outN.bad()) {
		cerr << "Failed to save '" << a2N.getGenome().getName() << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(outN.str() != out.str()) {
		cerr << "Loaded and re-saved GenomeAnno objects don't match" << endl;
		return EXIT_FAILURE;
	}

	/* GFF IO tests */
	/* write individually */
	out.str("");
	a1.write(out);
	a2.write(out);
	/* write as MetaGenome annotations */
	ostringstream outM;
	MetaGenomeAnno aM;
	aM.push_back(a1);
	aM.push_back(a2);
	aM.write(outM);
	if(out.str() != outM.str()) {
		cerr << "Individual and MetaGenome annotations don't match: -----" << endl
				<< out.str() << endl << "-----" << endl
				<< outM.str() << endl << "-----" << endl;
		return EXIT_FAILURE;
	}

	/* part 2, MetaGenomeAnno test */
	out.str("");
	MetaGenomeAnno mta1, mta2;
	mta1.push_back(a1);
	mta1.push_back(a2);
	mta2.push_back(a1N);
	mta2.push_back(a2N);

	MetaGenomeAnno mta = mta1 + mta2;

	mta.save(out);
	if(out.bad()) {
		cerr << "Failed to save MetaGenomeAnno: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	in.str(out.str());
	MetaGenomeAnno mtaN;
	mtaN.load(in);
	if(in.bad()) {
		cerr << "Failed to load MetaGenomeAnno: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	ostringstream gffO, gffON;
	mta.write(gffO);
	mtaN.write(gffON);

	if(gffON.str() != gffO.str()) {
		cerr << "Unmatched MetaGenone GFF annotation" << endl;
		cerr << gffO.str() << endl <<
				"-----" << endl <<
				gffON.str() << endl;
	}

	cout << "MetaGenomeAnno tests passed" << endl;
}
