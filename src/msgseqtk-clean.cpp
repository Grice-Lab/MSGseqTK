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
//============================================================================
// Name        : msgseqtk.cpp
// Author      : Qi Zheng
// Version     :
// Copyright   : GPL v3.0 Copyright (C) 2017  Qi Zheng
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
//#include <boost/random/mersenne_twister.hpp>
#include "MSGseqTK.h"
#include "EGUtil.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqTK;

/* program default values */
static const double DEFAULT_MIN_LOD = 0;
//static const double MAX_EVAL = 1;
static const int DEFAULT_STRAND = 3;
//static const double DEFAULT_INDEL_RATE = 0.01;
static const string ASSIGNMENT_BASIC_HEADER = "id\tdescription\tref_loglik\tbg_loglik\tLOD";
static const string ASSIGNMENT_SE_DETAIL = "ref_MEMS\tbg_MEMS";
static const string ASSIGNMENT_PE_DETAIL = "ref_MEMS_fwd\tref_MEMS_rev\tbg_MEMS_fwd\tbg_MEMS_rev";
static const int DEFAULT_NUM_THREADS = 1;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Clean/remove background (i.e. host contamination) reads from Metagenomics Shot-Gun (MSG) sequencing data,"
		 << " based on Maximal Exact Matched Seeds (MEMS) searches and Baysian inference" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif

	cerr << "Usage:    " << progName << "  <READ-FILE> [MATE-FILE] <-r|--ref REF-DB> <-b|--bg BG-DB> <-o READ-OUT> [-p MATE-OUT] [other-options]" << endl
		 << "READ-FILE  FILE                  : single-end/forward read file need to be cleaned" << ZLIB_SUPPORT << endl
		 << "MATE-FILE  FILE                  : mate/reverse read file need to be cleaned" << ZLIB_SUPPORT << endl
	     << "Options:    -r|--ref  STR        : name of reference/target database from which WGS reads need to be kept" << endl
		 << "            -b|--bg  STR         : name of background/host database from which WGS reads need to be removed" << endl
		 << "            -o  FILE             : output of cleaned single-end/forward reads" << ZLIB_SUPPORT << endl
		 << "            -p  FILE             : output of cleaned mate/reverse reads" << ZLIB_SUPPORT << endl
		 << "            -a  FILE             : write an additional TSV file with the detailed assignment information for each read" << endl
		 << "            --detail  FLAG       : write detailed information of MEMS in the assignment output, ignored if -a is not specified; this may lead to slower processing" << endl
		 << "            -L|--lod  DBL        : minimum log-odd required to determine a read/pair as reference vs. background [" << DEFAULT_MIN_LOD << "]" << endl
//		 << "            -e  DBL              : maximum e-value allowed to consider an MEM between datbase and a read as significance [" << MAX_EVAL << "]" << endl
		 << "            -s|--strand  INT     : read/pair strand to search, 1 for sense, 2 for anti-sense, 3 for both [" << DEFAULT_STRAND << "]" << endl
		 << "            -S|--seed  INT       : random seed used for determing MEMS, for debug only" << endl
#ifdef _OPENMP
		 << "            -p|--process INT     : number of threads/cpus used for parallel processing" << endl
#endif
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

/**
 * main function to process single-ended reads
 * @return 0 if success, return non-zero otherwise
 */
int main_SE(MetaGenome& refMtg, MetaGenome& bgMtg, FMIndex& refFmidx, FMIndex& bgFmidx,
		SeqIO& seqI, SeqIO& seqO, ofstream& assignOut, RNG& rng,
		bool withDetail, double minLod, int strand);

/**
 * main function to process paired-ended reads
 * @return 0 if success, return non-zero otherwise
 */
int main_PE(MetaGenome& refMtg, MetaGenome& bgMtg, FMIndex& refFmidx, FMIndex& bgFmidx,
		SeqIO& fwdI, SeqIO& revI, SeqIO& fwdO, SeqIO& revO, ofstream& assignOut, RNG& rng,
		bool withDetail, double minLod, int strand);

int main(int argc, char* argv[]) {
	/* variable declarations */
	string fwdInFn, revInFn;
	string refDB, bgDB;
	string fwdOutFn, revOutFn, assignFn;

	boost::iostreams::filtering_istream fwdIn, revIn;
	boost::iostreams::filtering_ostream fwdOut, revOut;

	ofstream assignOut;
	bool withDetail = false;

	double minLod = DEFAULT_MIN_LOD;
	int strand = DEFAULT_STRAND;
	unsigned seed = time(NULL); // using time as default seed
//	double maxIndelRate = DEFAULT_INDEL_RATE;
	int nThreads = DEFAULT_NUM_THREADS;

	typedef boost::random::mt11213b RNG; /* random number generator type */

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

	if(!(1 <= cmdOpts.numMainOpts() && cmdOpts.numMainOpts() <= 2)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	fwdInFn = cmdOpts.getMainOpt(0);
	if(cmdOpts.numMainOpts() == 2)
		revInFn = cmdOpts.getMainOpt(1);

	if(cmdOpts.hasOpt("-r"))
		refDB = cmdOpts.getOpt("-r");
	if(cmdOpts.hasOpt("--ref"))
		refDB = cmdOpts.getOpt("--ref");

	if(cmdOpts.hasOpt("-b"))
		bgDB = cmdOpts.getOpt("-b");
	if(cmdOpts.hasOpt("--bg"))
		bgDB = cmdOpts.getOpt("--bg");

	if(cmdOpts.hasOpt("-o"))
		fwdOutFn = cmdOpts.getOpt("-o");
	if(cmdOpts.hasOpt("-p"))
		revOutFn = cmdOpts.getOpt("-p");

	if(cmdOpts.hasOpt("-a"))
		assignFn = cmdOpts.getOpt("-a");

	if(cmdOpts.hasOpt("--detail"))
		withDetail = true;

	if(cmdOpts.hasOpt("-q"))
		minLod = ::atof(cmdOpts.getOptStr("-q"));

	if(cmdOpts.hasOpt("-s"))
		strand = ::atoi(cmdOpts.getOptStr("-s"));
	if(cmdOpts.hasOpt("--strand"))
		strand = ::atoi(cmdOpts.getOptStr("--strand"));

	if(cmdOpts.hasOpt("-S"))
		seed = ::atoi(cmdOpts.getOptStr("-S"));
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOptStr("--seed"));
	srand(seed); /* set seed */

#ifdef _OPENMP
	if(cmdOpts.hasOpt("-p"))
		nThreads = ::atoi(cmdOpts.getOptStr("-p"));
	if(cmdOpts.hasOpt("--process"))
		nThreads = ::atoi(cmdOpts.getOptStr("--process"));
#endif

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(refDB.empty()) {
		cerr << "reference DB must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(bgDB.empty()) {
		cerr << "background must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(fwdOutFn.empty()) {
		cerr << "-o must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(!revInFn.empty() && revOutFn.empty()) {
		cerr << "-p must be specified when reverse read file provided" << endl;
		return EXIT_FAILURE;
	}

	if(!(minLod >= 0)) {
		cerr << "-q must be non-negative" << endl;
		return EXIT_FAILURE;
	}

#ifdef _OPENMP
	if(!(nThreads > 0)) {
		cerr << "-p|--process must be positive" << endl;
		return EXIT_FAILURE;
	}
	omp_set_num_threads(nThreads);
#endif

	/* guess input seq format */
	string fmt = SeqIO::guessFormat(fwdInFn);
	if(!(fmt == "fasta" || fmt == "fastq")) {
		cerr << "Unrecognized sequence format for file '" << fwdInFn << "'" << endl;
		return EXIT_FAILURE;
	}

	bool isPaired = !revInFn.empty();

	/* open inputs */
#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(fwdInFn, GZIP_FILE_SUFFIX))
		fwdIn.push(boost::iostreams::gzip_decompressor());
	else if(StringUtils::endsWith(fwdInFn, BZIP2_FILE_SUFFIX))
		fwdIn.push(boost::iostreams::bzip2_decompressor());
	else { }
#endif
	fwdIn.push(boost::iostreams::file_source(fwdInFn));
	if(fwdIn.bad()) {
		cerr << "Unable to open '" << fwdInFn << "' " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(isPaired) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(revInFn, GZIP_FILE_SUFFIX))
			revIn.push(boost::iostreams::gzip_decompressor());
		else if(StringUtils::endsWith(revInFn, BZIP2_FILE_SUFFIX))
			revIn.push(boost::iostreams::bzip2_decompressor());
		else { }
#endif
		revIn.push(boost::iostreams::file_source(revInFn));
		if(revIn.bad()) {
			cerr << "Unable to open '" << revInFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* open SeqIO input */
	SeqIO fwdI(dynamic_cast<istream*>(&fwdIn), fmt);
	SeqIO revI(dynamic_cast<istream*>(&revIn), fmt);
	string refMtgFn = refDB + METAGENOME_FILE_SUFFIX;
	string refFmidxFn = refDB + FMINDEX_FILE_SUFFIX;
	string bgMtgFn = bgDB + METAGENOME_FILE_SUFFIX;
	string bgFmidxFn = bgDB + FMINDEX_FILE_SUFFIX;

	ifstream refMtgIn, bgMtgIn;
	ifstream refFmidxIn, bgFmidxIn;

	MetaGenome refMtg, bgMtg;
	FMIndex refFmidx, bgFmidx;

	refMtgIn.open(refMtgFn.c_str(), ios_base::binary);
	if(!refMtgIn.is_open()) {
		cerr << "Unable to open '" << refMtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	refFmidxIn.open(refFmidxFn.c_str(), ios_base::binary);
	if(!refFmidxIn.is_open()) {
		cerr << "Unable to open '" << refFmidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	bgMtgIn.open(bgMtgFn.c_str(), ios_base::binary);
	if(!bgMtgIn.is_open()) {
		cerr << "Unable to open '" << bgMtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	bgFmidxIn.open(bgFmidxFn.c_str(), ios_base::binary);
	if(!bgFmidxIn.is_open()) {
		cerr << "Unable to open '" << bgFmidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(fwdOutFn, GZIP_FILE_SUFFIX)) /* empty outFn won't match */
		fwdOut.push(boost::iostreams::gzip_compressor());
	else if(StringUtils::endsWith(fwdOutFn, BZIP2_FILE_SUFFIX)) /* empty outFn won't match */
		fwdOut.push(boost::iostreams::bzip2_compressor());
	else { }
#endif
	fwdOut.push(boost::iostreams::file_sink(fwdOutFn));
	if(fwdOut.bad()) {
		cerr << "Unable to write to '" << fwdOutFn << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(isPaired) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(revOutFn, GZIP_FILE_SUFFIX)) /* empty outFn won't match */
			revOut.push(boost::iostreams::gzip_compressor());
		else if(StringUtils::endsWith(revOutFn, BZIP2_FILE_SUFFIX)) /* empty outFn won't match */
			revOut.push(boost::iostreams::bzip2_compressor());
		else { }
#endif
		revOut.push(boost::iostreams::file_sink(revOutFn));
		if(revOut.bad()) {
			cerr << "Unable to write to '" << revOutFn << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!assignFn.empty()) {
		assignOut.open(assignFn.c_str());
		if(!assignOut.is_open()) {
			cerr << "Unable to write to '" << assignFn << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* open SeqIO output */
	SeqIO fwdO(dynamic_cast<ostream*>(&fwdOut), fmt);
	SeqIO revO(dynamic_cast<ostream*>(&revOut), fmt);

	/* load data */
	infoLog << "Loading reference MetaGenome info ..." << endl;
	loadProgInfo(refMtgIn);
	if(!refMtgIn.bad())
		refMtg.load(refMtgIn);
	if(refMtgIn.bad()) {
		cerr << "Unable to load reference MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Loading refrence FM-index ..." << endl;
	loadProgInfo(refFmidxIn);
	if(!refFmidxIn.bad())
		refFmidx.load(refFmidxIn);
	if(refFmidxIn.bad()) {
		cerr << "Unable to load reference FM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Loading background MetaGenome info ..." << endl;
	loadProgInfo(bgMtgIn);
	if(!bgMtgIn.bad())
		bgMtg.load(bgMtgIn);
	if(bgMtgIn.bad()) {
		cerr << "Unable to load background MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Loading background FM-index ..." << endl;
	loadProgInfo(bgFmidxIn);
	if(!bgFmidxIn.bad())
		bgFmidx.load(bgFmidxIn);
	if(bgFmidxIn.bad()) {
		cerr << "Unable to load background FM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* initiate RNG */
	RNG rng(seed);

	if(assignOut.is_open()) {
		writeProgInfo(assignOut, string(" read assignment generated by ") + argv[0]);
		assignOut << "# command: "<< cmdOpts.getCmdStr() << endl;
		string header = ASSIGNMENT_BASIC_HEADER;
		if(withDetail)
			header += "\t" + (!isPaired ? ASSIGNMENT_SE_DETAIL : ASSIGNMENT_PE_DETAIL);
		assignOut << header << endl;
	}

	/* main processing */
	if(!isPaired)
		return main_SE(refMtg, bgMtg, refFmidx, bgFmidx, fwdI, fwdO, assignOut, rng, withDetail, minLod, strand);
	else
		return main_PE(refMtg, bgMtg, refFmidx, bgFmidx, fwdI, revI, fwdO, revO, assignOut, rng, withDetail, minLod, strand);
}

int main_SE(MetaGenome& refMtg, MetaGenome& bgMtg, FMIndex& refFmidx, FMIndex& bgFmidx,
		SeqIO& seqI, SeqIO& seqO, ofstream& assignOut, RNG& rng,
		bool withDetail, double minLod, int strand) {
	infoLog << "Filtering input reads" << endl;
	/* search MEMS for each read */
#pragma omp parallel
	{
#pragma omp single
		{
			while(seqI.hasNext()) {
				PrimarySeq read = seqI.nextSeq();
				const string& id = read.getName();
				const string& desc = read.getDesc();
#pragma omp task
				{
					MEMS refMems = MEMS::sampleMEMS(&read, &refFmidx, rng, strand, withDetail);
					MEMS bgMems = MEMS::sampleMEMS(&read, &bgFmidx, rng, strand, withDetail);
					double refLoglik = refMems.loglik();
					double bgLoglik = bgMems.loglik();
					double lod = - refLoglik + bgLoglik;
					if(lod >= minLod)
#pragma omp critical(writeSeq)
						seqO.writeSeq(read);
					if(assignOut.is_open()) {
#pragma omp critical(writeAssign)
						{
							assignOut << id << "\t" << desc << "\t" << refLoglik << "\t" << bgLoglik << "\t" << lod;
							if(withDetail)
								assignOut << "\t" << refMems << "\t" << bgMems;
							assignOut << endl;
						}
					}
				} /* end task */
			} /* end each read */
		} /* end single */
#pragma omp taskwait
	} /* end parallel */
	return EXIT_SUCCESS;
}

int main_PE(MetaGenome& refMtg, MetaGenome& bgMtg, FMIndex& refFmidx, FMIndex& bgFmidx,
		SeqIO& fwdI, SeqIO& revI, SeqIO& fwdO, SeqIO& revO, ofstream& assignOut, RNG& rng,
		bool withDetail, double minLod, int strand) {
	infoLog << "Filtering input paired-end reads" << endl;
	/* search MEMS for each pair */
#pragma omp parallel
	{
#pragma omp single
		{
			while(fwdI.hasNext() && revI.hasNext()) {
				PrimarySeq fwdRead = fwdI.nextSeq();
				PrimarySeq revRead = revI.nextSeq();
				assert(fwdRead.getName() == revRead.getName());
				const string& id = fwdRead.getName();
				const string& desc = fwdRead.getDesc();
#pragma omp task
				{
					MEMS_PE refMems_pe = MEMS::sampleMEMS(&fwdRead, &revRead, &refFmidx, rng, strand, withDetail);
					MEMS_PE bgMems_pe = MEMS::sampleMEMS(&fwdRead, &revRead, &bgFmidx, rng, strand, withDetail);
					double refLoglik = MEMS::loglik(refMems_pe);
					double bgLoglik = MEMS::loglik(bgMems_pe);
					double lod = - refLoglik + bgLoglik;
					if(lod >= minLod) {
#pragma omp critical(writeSeq)
						{
							fwdO.writeSeq(fwdRead);
							revO.writeSeq(revRead);
						}
					}
					if(assignOut.is_open()) {
#pragma omp critical(writeAssign)
						{
							assignOut << id << "\t" << desc << "\t" << refLoglik << "\t" << bgLoglik << "\t" << lod;
							if(withDetail)
								assignOut << refMems_pe.first << "\t" << refMems_pe.second << "\t"
								<< bgMems_pe.first << "\t" << bgMems_pe.second;
							assignOut << endl;
						}
					}
				} /* end task */
			} /* end each pair */
		} /* end single */
#pragma omp taskwait
	} /* end parallel */
	return EXIT_SUCCESS;
}
