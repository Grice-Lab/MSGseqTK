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
#include <ctime>
#include <cmath>
#include <boost/algorithm/string.hpp> /* for boost string algorithms */
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include "EGUtil.h"
#include "MSGseqTK.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;
using namespace Eigen;

/** default program values */
static const int DEFAULT_READ_LEN = 100;
static const double DEFAULT_MEAN_INSERT = 300;
static const double DEFAULT_SD_INSERT = 50;
static const double DEFAULT_MIN_INSERT = 100;
static const double DEFAULT_MAX_INSERT = 500;
static const string READ_FMT = "fastq";
static const string READ_PREFIX = "r";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Generate simulated MetaGenomics Shot-gun sequencing (MGS) reads from a MSGseqTK database" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif

	cerr << "Usage:    " << progName << "  <DB-NAME> <-N NUM-READ> <-o READ-OUTFILE> [options]" << endl
		 << "DB-NAME  STR                     : MSGseqTK database name" << endl
		 << "Options:    -N  LONG             : # of simulated read" << endl
		 << "            -L  INT              : WGS read length [" << DEFAULT_READ_LEN << "]" << endl
		 << "            -o  FILE             : FASTQ output of simulated (forward) read" << ZLIB_SUPPORT << endl
		 << "            -p  FILE             : FASTQ output of simulated (reverse) read" << ZLIB_SUPPORT << endl
		 << "            -m|--mean  DBL       : mean length of WGS inserts [" << DEFAULT_MEAN_INSERT << "]" << endl
		 << "            -s|--sd  DBL         : sd of WGS inserts [" << DEFAULT_SD_INSERT << "]" << endl
		 << "            -M|--min  DBL        : min length of WGS inserts [" << DEFAULT_MIN_INSERT << "]" << endl
		 << "            -N|--max  DBL        : min length of WGS inserts [" << DEFAULT_MAX_INSERT << "]" << endl
		 << "            -S|--seed  INT       : random seed used for simulation, for debug purpose" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName;
	string fwdFn, revFn;
	ifstream mtgIn, fmidxIn;
	boost::iostreams::filtering_ostream fwdOut, revOut;
	SeqIO fwdSeqO, revSeqO;

	long N = 0;
	int readLen = DEFAULT_READ_LEN;
	double meanInsert = DEFAULT_MEAN_INSERT;
	double sdInsert = DEFAULT_SD_INSERT;
	double minInsert = DEFAULT_MIN_INSERT;
	double maxInsert = DEFAULT_MAX_INSERT;

	unsigned seed = time(NULL); // using time as default seed

	typedef boost::random::mt11213b RNG; /* random number generator type */
	typedef boost::random::uniform_int_distribution<size_t> LocDistrib; /* location distribution in metagenome/FM-index */
	typedef boost::random::normal_distribution<> SizeDistrib; /* read size distribution */
	//typedef boost::random::uniform_01<> GapDistrib; /* gap observing distribution */
//	typedef boost::random::discrete_distribution<int8_t> BaseDistrib; /* base (nucleotide) distribution */
//	typedef BaseDistrib::param_type BaseParam; /* base distribution parameters */

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

	if(cmdOpts.hasOpt("-N"))
		N = ::atol(cmdOpts.getOptStr("-N"));

	if(cmdOpts.hasOpt("-L"))
		readLen = ::atoi(cmdOpts.getOptStr("-L"));

	if(cmdOpts.hasOpt("-o"))
		fwdFn = cmdOpts.getOpt("-o");
	if(cmdOpts.hasOpt("-p"))
		revFn = cmdOpts.getOpt("-p");

	if(cmdOpts.hasOpt("-m"))
		meanInsert = ::atof(cmdOpts.getOptStr("-m"));
	if(cmdOpts.hasOpt("--mean"))
		meanInsert = ::atof(cmdOpts.getOptStr("--mean"));

	if(cmdOpts.hasOpt("-s"))
		sdInsert = ::atof(cmdOpts.getOptStr("-s"));
	if(cmdOpts.hasOpt("--sd"))
		sdInsert = ::atof(cmdOpts.getOptStr("--sd"));

	if(cmdOpts.hasOpt("-M"))
		minInsert = ::atof(cmdOpts.getOptStr("-M"));
	if(cmdOpts.hasOpt("--min"))
		minInsert = ::atof(cmdOpts.getOptStr("--min"));

	if(cmdOpts.hasOpt("-N"))
		maxInsert = ::atof(cmdOpts.getOptStr("-N"));
	if(cmdOpts.hasOpt("--max"))
		maxInsert = ::atof(cmdOpts.getOptStr("--max"));

	if(cmdOpts.hasOpt("-S"))
		seed = ::atoi(cmdOpts.getOptStr("-S"));
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOptStr("--seed"));

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(!(N > 0)) {
		cerr << "-N must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(readLen > 0)) {
		cerr << "-L must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(fwdFn.empty()) {
		cerr << "-o must be specificied" << endl;
		return EXIT_FAILURE;
	}
	if(meanInsert <= 0) {
		cerr << "-m|--mean must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(sdInsert < 0) {
		cerr << "-s|--sd must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(minInsert < 0) {
		cerr << "-M|--min must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(!(maxInsert >= minInsert)) {
		cerr << "-N|--max must be no less than -M|--min" << endl;
		return EXIT_FAILURE;
	}

	/* set dbName */
	string mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	string rfmFn = dbName + FMINDEX_FILE_SUFFIX;

	/* open inputs */
	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	fmidxIn.open(rfmFn.c_str(), ios_base::binary);
	if(!fmidxIn.is_open()) {
		cerr << "Unable to open '" << rfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(fwdFn, GZIP_FILE_SUFFIX))
		fwdOut.push(boost::iostreams::gzip_compressor());
	else if(StringUtils::endsWith(fwdFn, BZIP2_FILE_SUFFIX)) /* empty outFn won't match */
		fwdOut.push(boost::iostreams::bzip2_compressor());
	else { }
#endif

	fwdOut.push(boost::iostreams::file_sink(fwdFn));

	if(fwdOut.bad()) {
		cerr << "Unable to write to '" << fwdFn << "'" << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(!revFn.empty()) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(revFn, GZIP_FILE_SUFFIX))
			revOut.push(boost::iostreams::gzip_compressor());
		else if(StringUtils::endsWith(revFn, BZIP2_FILE_SUFFIX)) /* empty outFn won't match */
			revOut.push(boost::iostreams::bzip2_compressor());
		else { }
#endif

		revOut.push(boost::iostreams::file_sink(revFn));
		if(revOut.bad()) {
			cerr << "Unable to write to '" << revFn << "'" << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* prepare SeqIO */
	fwdSeqO.reset(dynamic_cast<ostream*>(&fwdOut), READ_FMT);
	if(!revFn.empty())
		revSeqO.reset(dynamic_cast<ostream*>(&revOut), READ_FMT);

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

	infoLog << "Reconstrucing metagenome sequence ..." << endl;
	const DNAseq& mtgSeq = fmidx.getSeq();
	const size_t L = mtgSeq.length();

	/* init random distributions */
	RNG rng(seed);
	LocDistrib loc_dist(0, L - 1);
	SizeDistrib size_dist(meanInsert, sdInsert);

	for(long n = 1; n <= N;) {
		/* sample start */
		long start, end; /* 0-based */
		start = loc_dist(rng);

		int len = static_cast<int>(::round(size_dist(rng)));
		if(len < minInsert)
			len = minInsert;
		if(len > maxInsert)
			len = maxInsert;
		end = start + len - 1;

		size_t startGenomeIdx = mtg.getGenomeIndex(start);
		size_t endGenomeIdx = mtg.getGenomeIndex(end);
		size_t startChromIdx = mtg.getChromIndex(start);
		size_t endChromIdx = mtg.getChromIndex(end);
		/* check whether start and end on the same genome and chromosome */
		if(startGenomeIdx == endGenomeIdx == startChromIdx == endChromIdx) {
			DNAseq insertSeq = mtgSeq.substr(start, len);
			string rname = READ_PREFIX + boost::lexical_cast<string>(n);
			string desc = "start=" + boost::lexical_cast<string>(start) +
					";end=" + boost::lexical_cast<string>(end) +
					";genome=" + mtg.getGenome(startGenomeIdx).getName() +
					";chrom=" + mtg.getGenome(startGenomeIdx).getChrom(startChromIdx).name;

			/* write fwd read */
			fwdSeqO.writeSeq(PrimarySeq(insertSeq.substr(0, readLen), rname, desc));
			/* write rev read, if requested */
			if(!revFn.empty())
				revSeqO.writeSeq(PrimarySeq(insertSeq.revcom().substr(0, readLen), rname, desc));
			n++;
		}
	}
}
