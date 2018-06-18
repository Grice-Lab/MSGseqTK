/*
 * msgseqclean-build.cpp
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
#include "MSGseqClean.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::MSGseqClean;
using namespace Eigen;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Build a database from one or more genomic sequence files" << endl;
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
		 << "            -l  FILE             : tab-delimited list with 1st field sample-names and 2nd field genome filenames/paths" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	vector<string> inFns;
	map<string, string> genomeFn2Name;
	string dbName, listFn;
	ifstream listIn;
	ofstream mtgOut, fmidxOut;
	const string fmt = "fasta";

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

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(dbName.empty()) {
		cerr << "-n must be specified" << endl;
		return EXIT_FAILURE;
	}

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

	/* set dbName */
	string mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	string fmidxFn = dbName + FMINDEX_FILE_SUFFIX;

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

	MetaGenome mtg;
	FMIndex fmidx;

	/* process each file */
	for(const string& inFn : inFns) {
		string genomeName = genomeFn2Name.at(inFn);

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
		infoLog << "Merging into database ... ";
		fmidx = FMIndex(genomeSeq) + fmidx; /* always use ther fresh object as lhs */

		assert(mtg.getSize() == fmidx.length());
		infoLog << " done. Currrent # of genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;
	}

	infoLog << "MetaGenomics database build. Total # of genomes: " << mtg.numGenomes() << " size: " << mtg.getSize() << endl;

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
