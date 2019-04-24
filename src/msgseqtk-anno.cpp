/*
 * msgseqtk-anno.cpp
 *  Annotate a pre-built MSGseqTK database
 *  Created on: Dec 5, 2018
 *  Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <limits>
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

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Annotate a pre-built database with optional 3rd-party GFF files for the genomes in it" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed GTF/GFF3 files";
	#endif

	cerr << "Usage:    " << progName << "  <DB> [options]" << endl
		 << "DB          STR                  : database name/prefix" << endl
		 << "Options:    -o  FILE             : write database annotation to FILE instead of the default file [DB" << GFF::GFF3_SUFFIX << "]" << endl
		 << "            -l  FILE             : tab-delimited annotation list with first field unique genome IDs and last field external GFF annotation filenames, any middle fields are ignored" << ZLIB_SUPPORT << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	map<string, string> genomeId2GffFn;

	string dbName;
	string listFn, mtgFn, gffOutFn;

	ifstream listIn, mtgIn;

	ofstream gffOut;

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
		cerr << "Error: <DB-NAME> required" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	dbName = cmdOpts.getMainOpt(0);

	if(cmdOpts.hasOpt("-o"))
		gffOutFn = cmdOpts.getOpt("-o");
	else
		gffOutFn = MetaGenomeAnno::getDBAnnoFn(dbName);

	if(cmdOpts.hasOpt("-l"))
		listFn = cmdOpts.getOpt("-l");

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */

	/* open inputs */
	if(!listFn.empty()) {
		listIn.open(listFn.c_str());
		if(!listIn.is_open()) {
			cerr << "Unable to open '" << listFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "Reading in genome names from '" << listFn << "'" << endl;
		string line;
		while(std::getline(listIn, line)) {
			if(line.empty() || line.front() == '#')
				continue;
			vector<string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			if(fields.size() >= 2) {
				string id = fields.front();
				string fn = fields.back();
				if(genomeId2GffFn.count(id) > 0) {
					warningLog << "Non-unique genome ID " << id << " found in " << listFn << ", ignore" << endl;
					continue;
				}
				genomeId2GffFn[Genome::formatName(id)] = fn; // use formatted genome id
			}
		}
		listIn.close();
		if(!genomeId2GffFn.empty())
			infoLog << "Found " << genomeId2GffFn.size() << " user-provided GFF annotation files" << endl;
	}

	/* open outputs */
	gffOut.open(gffOutFn.c_str());
	if(!gffOut.is_open()) {
		cerr << "Unable to write to '" << gffOutFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* load MetaGenome */
	MetaGenome mtg;

	infoLog << "Loading MetaGenome in database '" << dbName << "'" << endl;
	mtgFn = dbName + METAGENOME_FILE_SUFFIX;

	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open database file '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	loadProgInfo(mtgIn);
	if(!mtgIn.bad())
		mtg.load(mtgIn); // no need to load the sequence for annotation
	if(mtgIn.bad()) {
		cerr << "Unable to load '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* construct basic MetaGenomeAnno */
	MetaGenomeAnno mta(mtg);

	/* process each genome GFF file, if exists */
	for(const Genome& genome : mtg.getGenomes()) {
		const string& genomeId = genome.getId();
		if(genomeId2GffFn.count(genomeId) == 0)
			continue;

		const string& gffFn = genomeId2GffFn.at(genomeId);

		infoLog << "Processing external GFF annotation in '" << gffFn << "' for " << genome.displayId() << endl;

		GFF::Version extVer = GFF::UNK; /* GFF version for this gffFn */

		/* open external GFF file, and guess GFF version */
		boost::iostreams::filtering_istream gffIn;
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(gffFn, GZIP_FILE_SUFFIX)) {
			gffIn.push(boost::iostreams::gzip_decompressor());
			extVer = GFF::guessVersion(StringUtils::removeEnd(gffFn, GZIP_FILE_SUFFIX));
		}
		else if(StringUtils::endsWith(gffFn, BZIP2_FILE_SUFFIX)) {
			gffIn.push(boost::iostreams::bzip2_decompressor());
			extVer = GFF::guessVersion(StringUtils::removeEnd(gffFn, BZIP2_FILE_SUFFIX));
		}
		else {
			extVer = GFF::guessVersion(gffFn);
		}
#endif
		boost::iostreams::file_source gffSrc(gffFn);
		if(!gffSrc.is_open()) {
			cerr << "Unable to open GFF file '" << gffFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		gffIn.push(boost::iostreams::file_source(gffFn));

		if(extVer == GFF::UNK)
			warningLog << "Unable to guess GFF version from filename '" << gffFn << "', trying reading its content" << endl;

		mta.read(gffIn, genome, extVer);
		if(gffIn.bad()) {
			cerr << "Failed to read GFF file '" << gffFn << "' :" << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	} /* end each genome */

	/* output */
	infoLog << "Writing database annotations ..." << endl;
	MetaGenomeAnno::writeGFFHeader(gffOut, dbName, MetaGenomeAnno::ANNO_VER);
	mta.write(gffOut, mtg);
	if(gffOut.bad()) {
		cerr << "Unable to write GFF annotation: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Database annotation written into '" << gffOutFn << "'" << endl;
}
