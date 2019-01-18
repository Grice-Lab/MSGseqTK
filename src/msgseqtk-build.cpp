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

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;

//static const int DEFAULT_NUM_THREADS = 1;
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
		 << "Options:    -n  STR              : database name" << endl
		 << "            -l  FILE             : tab-delimited genome list with 1st field unique genome IDs, 2nd filed genome names, 3nd field genomic sequence filenames; if provided, <SEQ-FILE> options are ignored" << ZLIB_SUPPORT << endl
		 << "            -r|--update  STR     : update database based on this old DB, it can be the same name as -n, which will overwrite the old database" << endl
		 << "            -b|--block  INT      : block size (in Mbp) for building FM-index, larget block size will lead to faster but more memory usage algorithm [" << DEFAULT_BLOCK_SIZE << "]" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	map<string, string> genomeId2Name;
	map<string, string> genomeId2Fn;

	string dbName, oldDBName;
	string listFn, mtgFn, fmidxFn;

	ifstream listIn, mtgIn, fmidxIn;

	ofstream mtgOut, fmidxOut;
	const SeqIO::FORMAT fmt = SeqIO::FASTA; // always use fasta format

	int blockSize = DEFAULT_BLOCK_SIZE;

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

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(dbName.empty()) {
		cerr << "-n must be specified" << endl;
		return EXIT_FAILURE;
	}
	if(dbName == oldDBName)
		warningLog << "Warning: old database '" << oldDBName << "' will be overwritten!" << endl;

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
				genomeId2Name[id] = name;
				genomeId2Fn[id] = fn;
			}
		}
		listIn.close();
		infoLog << "Found " << genomeId2Fn.size() << " user-provided genome information" << endl;
	}
	if(genomeId2Fn.empty()) {
		cerr << "At least one genome file must be provided" << endl;
		return EXIT_FAILURE;
	}

	/* open output */
	mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	fmidxFn = dbName + FMINDEX_FILE_SUFFIX;

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

	/* try to open existing DB */
	if(!oldDBName.empty()) { /* is an update */
		infoLog << "Loading old database '" << oldDBName << "'" << endl;
		/* set oldDBName */
		mtgFn = oldDBName + METAGENOME_FILE_SUFFIX;
		fmidxFn = oldDBName + FMINDEX_FILE_SUFFIX;

		mtgIn.open(mtgFn.c_str(), ios_base::binary);
		if(!mtgIn.is_open()) {
			cerr << "Unable to open old database file '" << mtgFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		fmidxIn.open(fmidxFn.c_str(), ios_base::binary);
		if(!fmidxIn.is_open()) {
			cerr << "Unable to open old database file '" << fmidxFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		loadProgInfo(mtgIn);
		if(!mtgIn.bad())
			mtg.load(mtgIn);
		if(mtgIn.bad()) {
			cerr << "Unable to load '" << mtgFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		loadProgInfo(fmidxIn);
		if(!fmidxIn.bad())
			fmidx.load(fmidxIn);
		if(fmidxIn.bad()) {
			cerr << "Unable to load '" << fmidxFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		mtgIn.close();
		fmidxIn.close();
	}

	/* read all genomic files */
	size_t nProcessed = 0;
	for(const std::map<string, string>::value_type& entry : genomeId2Fn) {
		const string& genomeId = entry.first;
		const string& genomeFn = entry.second;
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

		genomeIn.push(boost::iostreams::file_source(genomeFn));
		if(genomeIn.bad()) {
			cerr << "Unable to open genome seq file '" << genomeFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* read in genome sequences */
		Genome genome(genomeId, genomeName);
		DNAseq genomeSeq;
		infoLog << "Reading genome " << genome.displayId() << endl;
		SeqIO seqI(&genomeIn, fmt);
		while(seqI.hasNext()) {
			const PrimarySeq& chr = seqI.nextSeq();
			const string& chrName = chr.getName();
			const DNAseq& chrSeq = chr.getSeq();
			debugLog << "  adding " << chrName << " with length " << chrSeq.length() << endl;
			genome.addChrom(Genome::Chrom(chrName, chrSeq.length()));
			genomeSeq += chrSeq;
			genomeSeq.push_back(DNAalphabet::GAP_BASE);
		}

		/* add this genome */
		infoLog << "  genome " << genome.displayId() << " added into database" << endl;
		mtg.addGenome(genome, genomeSeq);
	}

	/* update final index */
	infoLog << "Updating MetaGenome indices" << endl;
	mtg.updateIndex();
	if(dbName == oldDBName && nProcessed == 0) {
		infoLog << "No new genomes found, database not modified, quit updating" << endl;
		return EXIT_SUCCESS;
	}

	infoLog << "Building FM-index incrementally" << endl;
	DNAseq blockSeq;
	blockSeq.reserve(blockSize * MBP_UNIT);
	size_t k = 0;
	size_t nBlock = 0;
	for(size_t i = 0; i < mtg.numGenomes(); ++i) {
		blockSeq += mtg.getGenomeSeq(i);
		nBlock++;
		bool isLast = i == mtg.numGenomes() - 1;

		/* process block, if large enough */
		if(blockSeq.length() >= blockSize * MBP_UNIT || isLast) { /* last genome or full block */
			if(!isLast)
				infoLog << "Adding " << nBlock << " genomes in block " << ++k << " into FM-index" << endl;
			else
				infoLog << "Adding " << nBlock << " genomes in block " << ++k << " into FM-index and building final sampled Suffix-Array" << endl;
			blockSeq.erase(blockSeq.length() - 1); // remove GAP_BASE terminal
			blockSeq.reverse(); // reverse genome sequences
			fmidx = FMIndex(blockSeq, isLast) + fmidx; /* always use freshly built FMIndex as lhs */
			blockSeq.clear();
			nBlock = 0;
			infoLog << "Currrent # of bases in FM-index: " << fmidx.length() << endl;
		}
	}
	assert(mtg.size() == fmidx.length());

	/* save output */
	saveProgInfo(mtgOut);
	mtg.save(mtgOut);
	if(mtgOut.bad()) {
		cerr << "Unable to save MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MetaGenome info saved to '" << mtgFn << "'" << endl;

	saveProgInfo(fmidxOut);
	fmidx.save(fmidxOut);
	if(fmidxOut.bad()) {
		cerr << "Unable to save FM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "FM-index saved to '" << fmidxFn << "'" << endl;

	if(!oldDBName.empty())
		infoLog << "Database updated. Newly added # of genomes: " << nProcessed << " new size: " << mtg.size() << endl;
	infoLog << "Database built. Total # of genomes: " << mtg.numGenomes() << " size: " << mtg.size() << endl;
}
