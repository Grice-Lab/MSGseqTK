/*
 * msgseqclean-merge.cpp
 *
 *  Created on: Jun 1, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include "EGUtil.h"
#include "MSGseqTK.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Merge two or more databases into a new database" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <DB1> <DB2> [DB3 ...] <-n DBNAME> [options]" << endl
		 << "D1, D2, ...    STR               : database names need to be merged" << endl
		 << "Options:    -n  STR              : new name for the merged database" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	vector<string> inDBNames;
	string dbName;
	ofstream mtgOut, fmidxOut;

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

	if(!(cmdOpts.numMainOpts() >= 2)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	inDBNames = cmdOpts.getMainOpt();

	if(cmdOpts.hasOpt("-n"))
		dbName = cmdOpts.getOpt("-n");

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(dbName.empty()) {
		cerr << "-n must be specified" << endl;
		return EXIT_FAILURE;
	}

	/* set dbName */
	string mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	string fmidxFn = dbName + FMINDEX_FILE_SUFFIX;

	/* open outputs */
	mtgOut.open(mtgFn.c_str(), ios_base::out | ios_base::binary);
	if(!mtgOut.is_open()) {
		cerr << "Unable to write to '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	fmidxOut.open(fmidxFn.c_str(), ios_base::out | ios_base::binary);
	if(!fmidxOut.is_open()) {
		cerr << "Unable to write to '" << fmidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	MetaGenome mtg;
	FMIndex fmidx;

	/* process each database */
	for(const string& inDB : inDBNames) {
		string mtgInfn = inDB + METAGENOME_FILE_SUFFIX;
		string fmidxInfn = inDB + FMINDEX_FILE_SUFFIX;
		/* open DB files */
		ifstream mtgIn(mtgInfn.c_str(), ios_base::binary);
		if(!mtgIn.is_open()) {
			cerr << "Unable to open '" << mtgInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		ifstream fmidxIn(fmidxInfn.c_str(), ios_base::binary);
		if(!fmidxIn.is_open()) {
			cerr << "Unable to open '" << fmidxInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* read in genome sequence, concatenated with Ns */
		infoLog << "Loading database '" << inDB << "'" << endl;
		MetaGenome mtgPart;
		FMIndex fmidxPart;
		loadProgInfo(mtgIn);
		if(!mtgIn.bad())
			mtgPart.load(mtgIn);
		if(mtgIn.bad()) {
			cerr << "Unable to load '" << mtgInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		loadProgInfo(fmidxIn);
		if(!fmidxIn.bad())
			fmidxPart.load(fmidxIn);
		if(fmidxIn.bad()) {
			cerr << "Unable to load '" << fmidxInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* incremental update */
		infoLog << "Merging into new database ... ";
		mtg = mtgPart + mtg;
		fmidx = fmidxPart + fmidx;
		assert(mtg.getSize() + mtg.numGenomes() == fmidx.length());
		infoLog << " done. currrent size: " << mtg.getSize() << endl;
	}

	infoLog << "MetaGenomics database merged. # of Genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;

	infoLog << "Saving database files ..." << endl;
	/* write database files, all with prepend program info */
	saveProgInfo(mtgOut);
	mtg.save(mtgOut);
	if(mtgOut.bad()) {
		cerr << "Unable to save MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MetaGenome info saved" << endl;

	saveProgInfo(fmidxOut);
	fmidx.save(fmidxOut);
	if(fmidxOut.bad()) {
		cerr << "Unable to save RFM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "RFM-index saved" << endl;
}
