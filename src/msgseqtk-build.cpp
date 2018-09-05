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
 *
 *  Created on: May 23, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cassert>
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
using namespace Eigen;

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
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif

	cerr << "Usage:    " << progName << "  <SEQ-FILE1> [SEQ-FILE2 SEQ-FILE3 ...] <-n DBNAME> [options]" << endl
		 << "SEQ-FILE  FILE                   : genome sequence file with one file per-genome in FASTA format" << ZLIB_SUPPORT << endl
		 << "Options:    -n  STR              : database name" << endl
		 << "            -l  FILE             : tab-delimited list with 1st field genome-names and 2nd field genomic file paths/names" << endl
		 << "            -r|--update  STR     : update database based on this old DB, it can be the same name as -n, which will overwrite the old database" << endl
		 << "            -f  FLAG             : during building/updating genomes already exist with the same names will be ignored, set this flag to force adding them" << endl
		 << "            -g|--gff3  FLAG      : write an additional metagenome annotation file in GFF3 format" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	vector<string> inFns;
	map<string, string> genomeFn2Name;
	string dbName, oldDBName;
	string listFn, gffFn, mtgFn, fmidxFn;

	ifstream listIn, mtgIn, fmidxIn;

	ofstream mtgOut, fmidxOut, gffOut;
	const string fmt = "fasta";

	bool isForce = false;

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

	if(!(cmdOpts.numMainOpts() > 0)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	for(int i = 0; i < cmdOpts.numMainOpts(); ++i) {
		string fn = cmdOpts.getMainOpt(i);
		inFns.push_back(fn);
		genomeFn2Name[fn] = fn; /* use filename as genome name by default */
	}

	if(cmdOpts.hasOpt("-n"))
		dbName = cmdOpts.getOpt("-n");

	if(cmdOpts.hasOpt("-l"))
		listFn = cmdOpts.getOpt("-l");

	if(cmdOpts.hasOpt("-g") || cmdOpts.hasOpt("--gff3"))
		gffFn = dbName + UCSC::GFF::GFF3_SUFFIX;

	if(cmdOpts.hasOpt("-r"))
		oldDBName = cmdOpts.getOpt("-r");
	if(cmdOpts.hasOpt("--update"))
		oldDBName = cmdOpts.getOpt("--update");

	if(cmdOpts.hasOpt("-f"))
		isForce = true;

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
		int nRead = 0;
		infoLog << "Read in genome names from " << listFn << endl;
		inFns.clear(); /* clear inFiles */
		string line;
		while(std::getline(listIn, line)) {
			if(line.front() == '#')
				continue;
			vector<string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			if(fields.size() >= 2) {
				string name = fields[0];
				string fn = fields[1];
				if(genomeFn2Name.count(fn)) { /* this is an input file */
					inFns.push_back(fn);
					genomeFn2Name[fn] = name; /* update the sample name */
					nRead++;
				}
			}
		}
		listIn.close();
		infoLog << nRead << " user-provided genome names read" << endl;
	}
	if(inFns.empty()) {
		cerr << "At least one valid genome file must be provided" << endl;
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
		fmidxIn.open(fmidxFn.c_str(), ios_base::binary);

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

	/* process each file */
	size_t nProcessed = 0;
	for(const string& inFn : inFns) {
		string genomeName = genomeFn2Name.at(inFn);
		if(!isForce && mtg.hasGenome(genomeName)) {
			warningLog << "Genome '" << genomeName << "' already exists in the database, ignore" << endl;
			continue;
		}

		/* open genome file */
		boost::iostreams::filtering_istream genomeIn;
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(inFn, GZIP_FILE_SUFFIX))
			genomeIn.push(boost::iostreams::gzip_decompressor());
		else if(StringUtils::endsWith(inFn, BZIP2_FILE_SUFFIX))
			genomeIn.push(boost::iostreams::bzip2_decompressor());
		else { }
#endif

		genomeIn.push(boost::iostreams::file_source(inFn));
		if(genomeIn.bad()) {
			cerr << "Unable to open genome seq file '" << inFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* read in genome sequence, concatenated with Ns */
		infoLog << "Reading genome '" << genomeName << "'" << endl;

		Genome genome(genomeName);
		SeqIO seqI(&genomeIn, fmt);
		DNAseq genomeSeq;
		while(seqI.hasNext()) {
			PrimarySeq chr = seqI.nextSeq().reverse(); /* alwasy use reversed sequence */
			string chrName = chr.getName();
			DNAseq chrSeq = chr.getSeq();
			genome.addChrom(chrName, chrSeq.length());

			genomeSeq += chrSeq;
			genomeSeq.push_back(DNAalphabet::N); /* add a null terminal after each chrom */
		}

		/* incremental update backward */
		mtg.push_front(genome);
		infoLog << "Adding into database ... ";
		fmidx = FMIndex(genomeSeq) + fmidx; /* always use ther fresh object as lhs */

		assert(mtg.getSize() == fmidx.length());
		infoLog << " done. Currrent # of genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;
		nProcessed++;
	}

	if(oldDBName.empty() || nProcessed > 0) { /* original db not modified */
		infoLog << "Building the final SA ..." << endl;
		fmidx.buildSA();
	}

	if(oldDBName.empty())
		infoLog << "MetaGenomics database build. Total # of genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;
	else
		infoLog << "MetaGenomics database updated. Newly added # of genomes: " << nProcessed << " total # of genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;

	if(dbName != oldDBName || nProcessed > 0) { /* building new or updated database */
		/* set db file names */
		mtgFn = dbName + METAGENOME_FILE_SUFFIX;
		fmidxFn = dbName + FMINDEX_FILE_SUFFIX;
		/* open output files */
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

		if(!gffFn.empty()) {
			gffOut.open(gffFn.c_str());
			if(!gffOut.is_open()) {
				cerr << "Unable to write to '" << gffFn << "': " << ::strerror(errno) << endl;
				return EXIT_FAILURE;
			}
		}

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

		if(gffOut.is_open()) {
			mtg.writeGFF(gffOut, UCSC::GFF::GFF3, progName);
			infoLog << "GFF3 annotation file written" << endl;
		}
	}
	else
		infoLog << "Database was not modified" << endl;
}
