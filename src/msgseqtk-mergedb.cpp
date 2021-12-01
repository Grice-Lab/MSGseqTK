/*******************************************************************************
 * This file is part of MSGseqTK, a Metagenomics Shot-Gun sequencing ToolKit
 * for ultra-fast and accurate MSG-seq cleaning, mapping and more,
 * based on space-efficient FMD-index on entire collection of meta-genomics sequences.
 * Copyright (C) 2018, 2019  Qi Zheng
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
 * msgseqtk-merge.cpp
 *
 *  Created on: Jun 1, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <cassert>
#include "EGUtil.h"
#include "MSGseqTK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;
using UCSC::GFF;

static const int DEFAULT_NUM_THREADS = 1;

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
		 << "            --keep-order  FLAG   : keep the order of DB list, by default starting from the smallest DB for best merging speed" << endl
		 << "            --sample-rate  INT   : sample rate for the Suffix-Array (SA), recommend use smaller value means faster location search but more RAM/disk usage when the reference genomes are highly redundant [" << FMDIndex::SA_SAMPLE_RATE << "]" << endl
#ifdef _OPENMP
		 << "            -p|--process INT     : number of threads/cpus for parallel SA building [" << DEFAULT_NUM_THREADS << "]" << endl
#endif
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

/* sort MSGseqTK database according to their sizes */
int sort_DB(vector<string>& inDBNames);

int main(int argc, char* argv[]) {
	/* variable declarations */
	vector<string> inDBNames;
	string dbName;
	bool keepOrder = false;
	ofstream mtgOut, mgsOut, fmdidxOut;

	int saSampleRate = FMDIndex::SA_SAMPLE_RATE;
	int nThreads = DEFAULT_NUM_THREADS;

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

	if(cmdOpts.hasOpt("--keep-order"))
		keepOrder = true;

	if(cmdOpts.hasOpt("--sample-rate"))
		saSampleRate = ::atoi(cmdOpts.getOptStr("--sample-rate"));

#ifdef _OPENMP
	if(cmdOpts.hasOpt("-p"))
		nThreads = ::atoi(cmdOpts.getOptStr("-p"));
	if(cmdOpts.hasOpt("--process"))
		nThreads = ::atoi(cmdOpts.getOptStr("--process"));
#endif

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(dbName.empty()) {
		cerr << "-n must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(std::find(inDBNames.begin(), inDBNames.end(), dbName) != inDBNames.end()) {
		cerr << "new DB name '" << dbName << " cannot be the same as any old DB name, abort merging" << endl;
		return EXIT_FAILURE;
	}

	if(!(saSampleRate > 0)) {
		cerr << "--sample-rate must be positive" << endl;
		return EXIT_FAILURE;
	}

#ifdef _OPENMP
	if(!(nThreads > 0)) {
		cerr << "-p|--process must be positive" << endl;
		return EXIT_FAILURE;
	}
	omp_set_num_threads(nThreads);
#endif

	/* set dbName */
	string mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	string mgsFn = dbName + METAGENOME_SEQ_FILE_SUFFIX;
	string fmdidxFn = dbName + FMDINDEX_FILE_SUFFIX;

	/* open outputs */
	mtgOut.open(mtgFn.c_str(), ios_base::out | ios_base::binary);
	if(!mtgOut.is_open()) {
		cerr << "Unable to write to '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	mgsOut.open(mgsFn.c_str(), ios_base::out | ios_base::binary);
	if(!mgsOut.is_open()) {
		cerr << "Unable to write to '" << mgsFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	fmdidxOut.open(fmdidxFn.c_str(), ios_base::out | ios_base::binary);
	if(!fmdidxOut.is_open()) {
		cerr << "Unable to write to '" << fmdidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	MetaGenome mtg;
	FMDIndex fmdidx;

	if(!keepOrder) { /* sort DB names by their size */
		int flag = sort_DB(inDBNames);
		if(flag != 0) {
			cerr << "Failed to sort input databses, aborting" << endl;
			return EXIT_FAILURE;
		}
	}

	for(vector<string>::const_iterator inDB = inDBNames.begin(); inDB != inDBNames.end(); ++inDB) {
		bool isLast = inDB == inDBNames.end() - 1;
		string mtgInfn = *inDB + METAGENOME_FILE_SUFFIX;
		string mgsInfn = *inDB + METAGENOME_SEQ_FILE_SUFFIX;
		string fmdidxInfn = *inDB + FMDINDEX_FILE_SUFFIX;
		string gffInfn = *inDB + GFF::GFF3_SUFFIX;
		/* open DB files */
		ifstream mtgIn(mtgInfn.c_str(), ios_base::binary);
		if(!mtgIn.is_open()) {
			cerr << "Unable to open '" << mtgInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		ifstream mgsIn(mgsInfn.c_str(), ios_base::binary);
		if(!mgsIn.is_open()) {
			cerr << "Unable to open '" << mgsInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		ifstream fmdidxIn(fmdidxInfn.c_str(), ios_base::binary);
		if(!fmdidxIn.is_open()) {
			cerr << "Unable to open '" << fmdidxInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* read in database files */
		infoLog << "Loading database '" << *inDB << "'" << endl;
		MetaGenome mtgPart;
		FMDIndex fmdidxPart;
		loadProgInfo(mtgIn);
		if(!mtgIn.bad())
			mtgPart.load(mtgIn);
		if(mtgIn.bad()) {
			cerr << "Unable to load '" << mtgInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		mtgPart.loadSeq(mgsIn);
		if(mgsIn.bad()) {
			cerr << "Unable to load '" << mgsInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		loadProgInfo(fmdidxIn);
		if(!fmdidxIn.bad())
			fmdidxPart.load(fmdidxIn);
		if(fmdidxIn.bad()) {
			cerr << "Unable to load '" << fmdidxInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* incremental update */
		if(!isLast)
			fmdidxPart.clearSA(); /* do not store SA info for intermediate parts */

		infoLog << "Merging into new database" << endl;

		mtg += mtgPart;
		/* check genome redundancy */
		size_t nGenome = mtg.removeRedundantGenomes();
		if(nGenome > 0) {
			errorLog << "Error: found " << nGenome << " redundant genome ids in metagenome database '"
					<< *inDB << "', aborting" << endl;
			return EXIT_FAILURE;
		}

		fmdidx += fmdidxPart;
		assert(mtg.BDSize() == fmdidx.length());
		infoLog << "Currrent # of genomes: " << mtg.numGenomes() << " # of bases: " << fmdidx.length() << endl;
	}

	/* check metagenome-level chromosome redundancy again */
	size_t nChr = mtg.renameRedundantChroms();
	if(nChr > 0)
		warningLog << "Warning: renamed " << nChr << " chromosome names at metagenome level to avoid redundancy"
		<< endl << " you may need to modify your external GFF annotation files before using msgseqtk-anno tool" << endl;
	infoLog << "Updating MetaGenome index" << endl;
	mtg.updateIndex();

	infoLog << "Building final Suffix-Array" << endl;
	fmdidx.buildSA(saSampleRate);

	infoLog << "Database merged. # of genomes: " << mtg.numGenomes() << " # of bases: " << fmdidx.length() << endl;
	infoLog << "Saving database files ..." << endl;
	/* write database files, all with prepend program info */
	saveProgInfo(mtgOut);
	mtg.save(mtgOut);
	if(mtgOut.bad()) {
		cerr << "Unable to save MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MetaGenome info saved" << endl;
	mtg.saveSeq(mgsOut);
	if(mgsOut.bad()) {
		cerr << "Unable to save MetaGenome seq: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MetaGenome seq saved" << endl;

	saveProgInfo(fmdidxOut);
	fmdidx.save(fmdidxOut);
	if(fmdidxOut.bad()) {
		cerr << "Unable to save FMD-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "FMD-index saved" << endl;
}

int sort_DB(vector<string>& inDBNames) {
	infoLog << "Determing input DB size" << endl;
	unordered_map<string, int64_t> inDBSizes;
	for(const string& dbName : inDBNames) {
		string mtgInfn = dbName + METAGENOME_FILE_SUFFIX;
		/* open DB mtg file */
		ifstream mtgIn(mtgInfn.c_str(), ios_base::binary);
		if(!mtgIn.is_open()) {
			cerr << "Unable to open '" << mtgInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		/* read in mtg file */
		MetaGenome mtg;
		loadProgInfo(mtgIn);
		if(!mtgIn.bad())
			mtg.load(mtgIn);
		if(mtgIn.bad()) {
			cerr << "Unable to load '" << mtgInfn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* record dbSizes */
		int64_t dbSize = mtg.size();
		debugLog << "  Database " << dbName << " size: " << dbSize << endl;
		inDBSizes[dbName] = dbSize;
	}

	/* sort DB by sizes */
	std::sort(inDBNames.begin(), inDBNames.end(),
			[&] (const string& lhs, const string& rhs)->bool { return inDBSizes[lhs] < inDBSizes[rhs]; }
	);

	return 0;
}
