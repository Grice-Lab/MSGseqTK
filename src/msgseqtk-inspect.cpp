/*******************************************************************************
 * This file is part of MSGseqTK, a Metagenomics Shot-Gun sequencing ToolKit
 * for ultra-fast and accurate MSG-seq cleaning, mapping and more,
 * based on space-efficient FM-index on entire collection of meta-genomics sequences.
 * Copyright (C) 2018  Qi Zheng
 *
 * MSGseqTK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MSGseqTK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/*
 * msgseqtk-inspect.cpp
 *
 *  Created on: Jun 1, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "EGUtil.h"
#include "MSGseqTK.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Check and verify the content of a previously built database" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <DB-NAME> [options]" << endl
		 << "DB-NAME    STR                   : database name (prefix)" << endl
		 << "Options:    -l  FILE             : write the genome names included in this database to FILE" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName;
	string listFn;
	ifstream mtgIn, fmidxIn, gffIn;
	ofstream listOut;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.empty() || cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printIntro();
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.hasOpt("--version")) {
		printVersion(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.numMainOpts() != 1) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	dbName = cmdOpts.getMainOpt(0);

	if(cmdOpts.hasOpt("-l"))
		listFn = cmdOpts.getOpt("-l");

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	/* set dbName */
	string mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	string fmidxFn = dbName + FMINDEX_FILE_SUFFIX;
	string gffFn = dbName + GFF::GFF3_SUFFIX;

	/* open inputs */
	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	fmidxIn.open(fmidxFn.c_str(), ios_base::binary);
	if(!fmidxIn.is_open()) {
		cerr << "Unable to open '" << fmidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	gffIn.open(gffFn.c_str());
	if(!gffIn.is_open()) {
		cerr << "Unable to open '" << gffFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
	if(!listFn.empty()) {
		listOut.open(listFn.c_str());
		if(!listOut.is_open()) {
			cerr << "Unable to write to '" << listFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* load data */
	MetaGenome mtg;
	FMIndex fmidx;

	infoLog << "Loading MetaGenome info ..." << endl;
	loadProgInfo(mtgIn);
	if(!mtgIn.bad())
		mtg.load(mtgIn);
	if(mtgIn.bad()) {
		cerr << "Unable to load MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Loading FM-index ..." << endl;
	loadProgInfo(fmidxIn);
	if(!fmidxIn.bad())
		fmidx.load(fmidxIn);
	if(fmidxIn.bad()) {
		cerr << "Unable to load FM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	MetaGenomeAnno mtgAnno(mtg);
	infoLog << "Loading MetaGenome annotation ..." << endl;
	string db;
	GFF::Version ver;
	MetaGenomeAnno::readGFFHeader(gffIn, db, ver);
	if(!(db == dbName && ver == GFF::GFF3)) {
		cerr << "Content from database annotation file '" << gffFn << " doesn't match its name" << endl;
		return EXIT_FAILURE;
	}
	mtgAnno.read(gffIn);
	if(gffIn.bad()) {
		cerr << "Unable to read MetaGenome annotation: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	cout << "MetaGenome info: # of genomes: " << mtg.numGenomes() << " size: " << mtg.size() << endl;
	cout << "FM-index info: length: " << fmidx.length() << " baseCount: ";
	for(int8_t b = DNAalphabet::A; b < DNAalphabet::SIZE; ++b)
		cout << " " << DNAalphabet::decode(b) << ":" << fmidx.getBaseCount(b);
	cout << endl;
	cout << "# of annotated genomes: " << mtgAnno.numAnnotated() <<
			" # of total annotations: " << mtgAnno.numAnnotations() << endl;

	if(listOut.is_open()) {
		for(const Genome& genome : mtg.getGenomes())
			listOut << genome.getId() << "\t" << genome.getName() << endl;
	}
}
