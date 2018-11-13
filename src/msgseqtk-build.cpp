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
static const int DEFAULT_BLOCK_SIZE = 1000;
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
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed files";
	#endif

	cerr << "Usage:    " << progName << "  <SEQ-FILE1> [SEQ-FILE2 SEQ-FILE3 ...] <-n DBNAME> [options]" << endl
		 << "SEQ-FILE  FILE                   : genome sequence file with one file per-genome in FASTA format" << ZLIB_SUPPORT << endl
		 << "Options:    -n  STR              : database name" << endl
		 << "            -l  FILE             : tab-delimited genome list with 1st field unique genome IDs, 2nd filed genome names, 3nd field genomic sequence file paths, and an optional 4th field with genomic GFF annotation file paths; if provided, <SEQ-FILE> options are ignored" << ZLIB_SUPPORT << endl
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
	map<string, string> genomeId2GffFn;

	string dbName, oldDBName;
	string oldGFFRecords;
	string listFn, gffFn, mtgFn, fmidxFn;

	ifstream listIn, mtgIn, fmidxIn, gffIn;

	ofstream mtgOut, fmidxOut, gffOut;
	const string fmt = "fasta";

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

	gffFn = dbName + UCSC::GFF::GFF3_SUFFIX;

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
				if(fields.size() > 3)
					genomeId2GffFn[id] = fields[3];
			}
		}
		listIn.close();
		infoLog << "Found " << genomeId2Fn.size() << " user-provided genome information" << endl;
		if(!genomeId2GffFn.empty())
			infoLog << "Found " << genomeId2GffFn.size() << " user-provided GFF annotation files" << endl;
	}
	if(genomeId2Fn.empty()) {
		cerr << "At least one genome file must be provided" << endl;
		return EXIT_FAILURE;
	}

	MetaGenome mtg;
	FMIndex fmidx;
	MetaGenomeAnno mtgAnno;

	/* try to open existing DB */
	if(!oldDBName.empty()) { /* is an update */
		infoLog << "Loading old database '" << oldDBName << "'" << endl;
		/* set oldDBName */
		mtgFn = oldDBName + METAGENOME_FILE_SUFFIX;
		fmidxFn = oldDBName + FMINDEX_FILE_SUFFIX;
		gffFn = oldDBName + GFF::GFF3_SUFFIX;

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
		gffIn.open(gffFn.c_str());
		if(!gffIn.is_open()) {
			cerr << "Unable to open old database annotation file '" << gffFn << "': " << ::strerror(errno) << endl;
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

		string db;
		GFF::Version ver;
		MetaGenomeAnno::readGFFHeader(gffIn, db, ver);
		if(db == oldDBName && ver == GFF::GFF3) {
			/* copy old records */
			oldGFFRecords = MetaGenomeAnno::readAll(gffIn);
			debugLog << "old GFF records copied from '" << gffFn << "'" << endl;
		}
		else
			warningLog << "Content from old database annotation file '" << gffFn << " doesn't match its name, ignored" << endl;

		mtgIn.close();
		fmidxIn.close();
		gffIn.close();
	}

	/* process each file */
	size_t nProcessed = 0;
	DNAseq blockSeq;
	vector<Genome> blockGenomes;
	int k = 0;
	for(const std::pair<string, string>& entry : genomeId2Fn) {
		string genomeId = entry.first;
		string genomeFn = entry.second;
		string genomeName = genomeId2Name[genomeId];

		if(mtg.hasGenome(genomeId)) {
			warningLog << "Genome " << genomeId << " (" << genomeName << ") already exists in the database, ignore" << endl;
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

		/* read in genome sequence, concatenated with Ns */
		infoLog << "Reading genome " << genomeId << " (" << genomeName << ")" << endl;

		Genome genome(genomeId, genomeName);
		SeqIO seqI(&genomeIn, fmt);

		while(seqI.hasNext()) {
			const PrimarySeq& chr = seqI.nextSeq();
			string chrName = chr.getName();
			DNAseq chrSeq = chr.getSeq().reverse(); /* reverse seq at each chrom */
			debugLog << "  adding " << chrName << " with length " << chrSeq.length() << endl;
			genome.addChrom(chrName, chrSeq.length());
			if(!blockSeq.empty())
				blockSeq.push_back(DNAalphabet::N); /* add an N terminal */
			blockSeq += chrSeq; /* N terminated chromosomes */
		}

		blockGenomes.push_back(genome); /* add this genome to the block */

		/* process external GFF file, if exists */
		if(genomeId2GffFn.count(genomeId)) {
			GenomeAnno anno(genome);
			gffFn = genomeId2GffFn[genomeId];
			GFF::Version extVer = GFF::UNK; /* GFF version for this gffFn */
			infoLog << "  Reading external GFF annotation from '" << gffFn << "'" << endl;

			/* open external GFF file, and guess GFF version */
			boost::iostreams::filtering_istream gffIn;
#ifdef HAVE_LIBZ
			if(StringUtils::endsWith(gffFn, GZIP_FILE_SUFFIX)) {
				gffIn.push(boost::iostreams::gzip_decompressor());
				extVer = GFF::guessVersion(StringUtils::removeEnd(static_cast<const string&>(gffFn), GZIP_FILE_SUFFIX));
			}
			else if(StringUtils::endsWith(gffFn, BZIP2_FILE_SUFFIX)) {
				gffIn.push(boost::iostreams::bzip2_decompressor());
				extVer = GFF::guessVersion(StringUtils::removeEnd(static_cast<const string&>(gffFn), BZIP2_FILE_SUFFIX));
			}
			else {
				extVer = GFF::guessVersion(gffFn);
			}
#endif
			gffIn.push(boost::iostreams::file_source(gffFn));

			if(extVer == GFF::UNK)
				extVer = GFF::guessVersion(gffFn);

			if(extVer != GFF::UNK) {
				if(!gffIn.bad()) {
					anno.read(gffIn, extVer);
					debugLog << "  read in " << anno.numAnnotated() << " external annotations" << endl;
					mtgAnno.push_back(anno);
				}
				else
					warningLog << "Unable to open external GFF file '" << gffFn << "' " << ::strerror(errno) << ", ignore" << endl;
			}
			else
				warningLog << "Unable to determine the GFF version of file '" << gffFn << "', ignore" << endl;
		} /* end processing external GFF file */

		/* process block, if large enough */
		bool isLast = genomeId2Fn.upper_bound(genomeId) == genomeId2Fn.end();
		if(blockSeq.length() >= blockSize * MBP_UNIT || isLast) { /* last genome or full block */
			if(!isLast)
				infoLog << "Adding " << blockGenomes.size() << " genomes in block " << ++k << " into database" << endl;
			else
				infoLog << "Adding " << blockGenomes.size() << " genomes in block " << ++k << " into database and building final sampled Suffix-Array" << endl;
			mtg.prepend(blockGenomes);
			fmidx = FMIndex(blockSeq, isLast) + fmidx; /* always use freshly built FMIndex as lhs */
			assert(mtg.size() == fmidx.length());
			blockGenomes.clear();
			blockSeq.clear();
			infoLog << "Currrent # of genomes: " << mtg.numGenomes() << " # of bases: " << fmidx.length() << endl;
		}
		/* incremental update backward */
		nProcessed++;
	}

	if(!oldDBName.empty())
		infoLog << "MetaGenomics database updated. Newly added # of genomes: " << nProcessed << endl;
	infoLog << "MetaGenomics database build. Total # of genomes: " << mtg.numGenomes() << " size: " << mtg.size() << endl;

	if(dbName != oldDBName || nProcessed > 0) { /* building new or updated database */
		/* set db file names */
		mtgFn = dbName + METAGENOME_FILE_SUFFIX;
		fmidxFn = dbName + FMINDEX_FILE_SUFFIX;
		gffFn = dbName + UCSC::GFF::GFF3_SUFFIX;
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

		gffOut.open(gffFn.c_str());
		if(!gffOut.is_open()) {
			cerr << "Unable to write to '" << gffFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
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
			cerr << "Unable to save FM-index: " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "FM-index saved" << endl;

		/* write MetaGenome annotations, and insert existing annotation if exist */
		MetaGenomeAnno::writeGFFHeader(gffOut, dbName, GenomeAnno::FORMAT);
		mtgAnno.write(gffOut);
		gffOut << oldGFFRecords; /* append all old records */
		infoLog << "GFF3 annotation file written" << endl;
	}
	else
		infoLog << "Database was not modified" << endl;
}
