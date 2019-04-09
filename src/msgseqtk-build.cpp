/*******************************************************************************
 * This file is part of MSGseqTK, a Metagenomics Shot-Gun sequencing ToolKit
 * for ultra-fast and accurate MSG-seq cleaning, mapping and more,
 * based on space-efficient FMD-index on entire collection of meta-genomics sequences.
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
 * msgseqtk-build.cpp
 *  build a MSGseqTK database
 *  Created on: May 23, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <limits>
#include <boost/algorithm/string.hpp> /* for boost string algorithms */
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
#include "EGUtil.h"
#include "MSGseqTK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;

static const int DEFAULT_NUM_THREADS = 1;
static const int DEFAULT_BLOCK_SIZE = 2000;
static const size_t MBP_UNIT = 1000000;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Build/update a database from one or more genomic sequence files" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed sequence files";
	#endif

	cerr << "Usage:    " << progName << "  <SEQ-FILE1> [SEQ-FILE2 SEQ-FILE3 ...] <-n DBNAME> [options]" << endl
		 << "SEQ-FILE  FILE                   : genome sequence file with one file per-genome in FASTA format" << ZLIB_SUPPORT << endl
		 << "Options:    -n  STR              : database name/prefix" << endl
		 << "            -l  FILE             : tab-delimited genome list with 1st field unique genome IDs, 2nd filed genome names, 3nd field genomic sequence filenames; if provided, <SEQ-FILE> options are ignored" << ZLIB_SUPPORT << endl
		 << "            -r|--update  STR     : update database based on this old DB, it can be the same name as -n, which will overwrite the old database" << endl
		 << "            -b|--block  INT      : block size (in Mbp) for building FMD-index, larget block size will lead to faster but more memory usage algorithm [" << DEFAULT_BLOCK_SIZE << "]" << endl
#ifdef _OPENMP
		 << "            -p|--process INT     : number of threads/cpus for parallel processing, only used in building SA and merging BWTs [" << DEFAULT_NUM_THREADS << "]" << endl
#endif
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	vector<string> genomeIds; // genome ids in original order
	map<string, string> genomeId2Name;
	map<string, string> genomeId2Fn;

	string dbName, oldDBName;
	string listFn, mtgFn, fmdidxFn;

	ifstream listIn, mtgIn, fmdidxIn;

	ofstream mtgOut, fmdidxOut;
	const SeqIO::FORMAT fmt = SeqIO::FASTA; // always use fasta format

	int blockSize = DEFAULT_BLOCK_SIZE;
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

	for(int i = 0; i < cmdOpts.numMainOpts(); ++i) {
		string fn = cmdOpts.getMainOpt(i);
		/* use filename as ID and Name */
		genomeIds.push_back(fn); // fn as id
		genomeId2Name[fn] = fn;
		genomeId2Fn[fn] = fn;
	}

	if(cmdOpts.hasOpt("-n"))
		dbName = cmdOpts.getOpt("-n");

	if(cmdOpts.hasOpt("-l"))
		listFn = cmdOpts.getOpt("-l");

	if(cmdOpts.hasOpt("-r"))
		oldDBName = cmdOpts.getOpt("-r");
	if(cmdOpts.hasOpt("--update"))
		oldDBName = cmdOpts.getOpt("--update");

	if(cmdOpts.hasOpt("-b"))
		blockSize = ::atoi(cmdOpts.getOptStr("-b"));
	if(cmdOpts.hasOpt("--block"))
		blockSize = ::atoi(cmdOpts.getOptStr("--block"));

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
	if(dbName == oldDBName)
		warningLog << "Warning: old database '" << oldDBName << "' will be overwritten!" << endl;

#ifdef _OPENMP
	if(!(nThreads > 0)) {
		cerr << "-p|--process must be positive" << endl;
		return EXIT_FAILURE;
	}
	omp_set_num_threads(nThreads);
#endif

	/* open inputs */
	if(!listFn.empty()) {
		listIn.open(listFn.c_str());
		if(!listIn.is_open()) {
			cerr << "Unable to open '" << listFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "Reading in genome names from '" << listFn << "'" << endl;
		genomeId2Name.clear();
		genomeId2Fn.clear();
		string line;
		while(std::getline(listIn, line)) {
			if(line.empty() || line.front() == '#')
				continue;
			vector<string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			if(fields.size() >= 3) {
				string id = fields[0];
				string name = fields[1];
				string fn = fields[2];
				if(genomeId2Name.count(id)) {
					warningLog << "Non-unique genome ID " << id << " found in " << listFn << ", ignore" << endl;
					continue;
				}
				genomeIds.push_back(id);
				genomeId2Name[id] = name;
				genomeId2Fn[id] = fn;
			}
		}
		listIn.close();
		infoLog << "Found " << genomeId2Fn.size() << " user-provided genome information" << endl;
	}
	if(genomeIds.empty()) {
		cerr << "At least one genome file must be provided" << endl;
		return EXIT_FAILURE;
	}

	/* open output */
	mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	fmdidxFn = dbName + FMDINDEX_FILE_SUFFIX;

	mtgOut.open(mtgFn.c_str(), ios_base::out | ios_base::binary);
	if(!mtgOut.is_open()) {
		cerr << "Unable to write to '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	fmdidxOut.open(fmdidxFn.c_str(), ios_base::out | ios_base::binary);
	if(!fmdidxOut.is_open()) {
		cerr << "Unable to write to '" << fmdidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	MetaGenome mtg;
	FMDIndex fmdidx;

	/* try to open existing DB */
	if(!oldDBName.empty()) { /* is an update */
		infoLog << "Loading old database '" << oldDBName << "'" << endl;
		/* set oldDBName */
		mtgFn = oldDBName + METAGENOME_FILE_SUFFIX;
		fmdidxFn = oldDBName + FMDINDEX_FILE_SUFFIX;

		mtgIn.open(mtgFn.c_str(), ios_base::binary);
		if(!mtgIn.is_open()) {
			cerr << "Unable to open old database file '" << mtgFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		fmdidxIn.open(fmdidxFn.c_str(), ios_base::binary);
		if(!fmdidxIn.is_open()) {
			cerr << "Unable to open old database file '" << fmdidxFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		loadProgInfo(mtgIn);
		if(!mtgIn.bad())
			mtg.load(mtgIn);
		if(mtgIn.bad()) {
			cerr << "Unable to load '" << mtgFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		loadProgInfo(fmdidxIn);
		if(!fmdidxIn.bad())
			fmdidx.load(fmdidxIn);
		if(fmdidxIn.bad()) {
			cerr << "Unable to load '" << fmdidxFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		mtgIn.close();
		fmdidxIn.close();
	}

	/* read all genomic files */
	size_t nProcessed = 0;
	for(const string& genomeId : genomeIds) {
		const string& genomeFn = genomeId2Fn[genomeId];
		const string& genomeName = genomeId2Name[genomeId];

		if(mtg.hasGenome(genomeId)) {
			warningLog << "Genome " << Genome::displayId(genomeId, genomeName) << " already exists in the database, ignore" << endl;
			continue;
		}

		/* open genome file */
		boost::iostreams::filtering_istream genomeIn;
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(genomeFn, GZIP_FILE_SUFFIX))
			genomeIn.push(boost::iostreams::gzip_decompressor());
		else if(StringUtils::endsWith(genomeFn, BZIP2_FILE_SUFFIX))
			genomeIn.push(boost::iostreams::bzip2_decompressor());
		else { }
#endif
		boost::iostreams::file_source genomeSrc(genomeFn);
		if(!genomeSrc.is_open()) {
			cerr << "Unable to open genome seq file '" << genomeFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		genomeIn.push(genomeSrc);

		/* read in genome sequences */
		Genome genome(genomeId, genomeName);
		infoLog << "Reading genome " << genome.displayId() << endl;
		SeqIO seqI(&genomeIn, fmt);
		while(seqI.hasNext()) {
			const PrimarySeq& chr = seqI.nextSeq();
			const string& chrName = chr.getName();
			const DNAseq& chrSeq = chr.getSeq();
			debugLog << "  adding " << chrName << " with length " << chrSeq.length() << endl;
			genome.addChrom(chrName, chrSeq);
		}

		/* add this genome */
		mtg.addGenome(genome);
		infoLog << "  genome " << genome.displayId() << " added" << endl;
	}

	/* update final index */
	infoLog << "Updating MetaGenome indices" << endl;
	mtg.updateIndex();
	if(dbName == oldDBName && nProcessed == 0) {
		infoLog << "No new genomes found, database not modified, quit updating" << endl;
		return EXIT_SUCCESS;
	}
	const size_t mtgSize = mtg.size();
	const size_t NG = mtg.numGenomes();

	/* save MetaGenome */
	saveProgInfo(mtgOut);
	mtg.save(mtgOut);
	if(mtgOut.bad()) {
		cerr << "Unable to save MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MetaGenome of " << mtgSize << " bases saved to '" << mtgFn << "'" << endl;

	infoLog << "Building FMD-index incrementally" << endl;
	DNAseq blockSeq;
	blockSeq.reserve(blockSize * MBP_UNIT);
	size_t k = 0;
	size_t nBlock = 0;
	while(!mtg.empty()) {
		blockSeq = dna::toBasic(mtg.topChrom().getBDSeq()) + blockSeq; // update seq
		mtg.popChrom(); // pop the last chrom

		nBlock++;
		bool isFirst = mtg.empty(); // flag whether the first genome

		/* process block, if large enough */
		if(blockSeq.length() >= blockSize * MBP_UNIT || isFirst) { /* first genome or full block */
			if(!isFirst)
				infoLog << "Adding " << nBlock << " chroms in block " << ++k << " into FMD-index" << endl;
			else
				infoLog << "Adding " << nBlock << " chroms in block " << ++k << " into FMD-index and building final sampled Suffix-Array" << endl;
			assert(blockSeq.back() == DNAalphabet::GAP_BASE);
			blockSeq.pop_back(); // remove GAP_BASE terminal
			fmdidx = FMDIndex(blockSeq, isFirst) + fmdidx; /* always use freshly built FMDIndex as lhs */
			blockSeq.clear();
			nBlock = 0;
			infoLog << "Currrent # of bases in FMD-index: " << fmdidx.length() << endl;
		}
	}
	assert(mtgSize == fmdidx.length());

	/* save FMDIndex */
	saveProgInfo(fmdidxOut);
	fmdidx.save(fmdidxOut);
	if(fmdidxOut.bad()) {
		cerr << "Unable to save FMD-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "FMD-index saved to '" << fmdidxFn << "'" << endl;

	if(!oldDBName.empty())
		infoLog << "Database updated. Newly added # of genomes: " << nProcessed << " new size: " << mtgSize << endl;
	infoLog << "Database built. Total # of genomes: " << NG << " size: " << mtgSize << endl;
}
