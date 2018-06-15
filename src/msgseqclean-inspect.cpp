/*
 * msgseqclean-inspect.cpp
 *
 *  Created on: Jun 1, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include "EGUtil.h"
#include "MSGseqClean.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqClean;

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
	ifstream mtgIn, rfmIn;
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
	string rfmFn = dbName + FMINDEX_FILE_SUFFIX;

	/* open inputs */
	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	rfmIn.open(rfmFn.c_str(), ios_base::binary);
	if(!rfmIn.is_open()) {
		cerr << "Unable to open '" << rfmFn << "': " << ::strerror(errno) << endl;
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
	loadProgInfo(rfmIn);
	if(!rfmIn.bad())
		fmidx.load(rfmIn);
	if(rfmIn.bad()) {
		cerr << "Unable to load FM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	cout << "MetaGenome info: # of genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;
	cout << "FM-index info: length: " << fmidx.length() << endl;
	for(int8_t b = DNAalphabet::A; b < DNAalphabet::SIZE; ++b)
		cout << " " << DNAalphabet::decode(b) << ":" << fmidx.getBaseCount(b);
	cout << endl;

	if(listOut.is_open()) {
		for(const string& genomeName : mtg.getGenomeNames())
			listOut << genomeName << endl;
	}
}
