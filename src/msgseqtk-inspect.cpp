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
 * msgseqtk-inspect.cpp
 *
 *  Created on: Jun 1, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <string>
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

static const string GENOME_TSV_HEADER = "genomeId\tgenomeName\tsize";
static const string CHROM_TSV_HEADER = "chromName\tsize\tgenomeId\tgenomeName";

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
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed sequence files";
	#endif

	cerr << "Usage:    " << progName << "  <DB-NAME> [options]" << endl
		 << "DB-NAME    STR                   : database name/prefix" << endl
		 << "Options:    -g|--genome  FILE    : write the genome names and length in this database in TSV fomrat" << endl
		 << "            -c|--chrom  FILE     : write the chromosome names and length in this database in TSV format" << endl
		 << "            -s|--seq  FILE       : write the genome sequences in this database " << ZLIB_SUPPORT << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName;
	string genomeFn, chromFn, seqFn, mtgFn, fmdidxFn;
	ifstream mtgIn, fmdidxIn;
	ofstream genomeOut, chromOut;
	boost::iostreams::filtering_ostream seqOut;

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

	if(cmdOpts.hasOpt("-g"))
		genomeFn = cmdOpts.getOpt("-g");
	if(cmdOpts.hasOpt("--genome"))
		genomeFn = cmdOpts.getOpt("--genome");

	if(cmdOpts.hasOpt("-c"))
		chromFn = cmdOpts.getOpt("-c");
	if(cmdOpts.hasOpt("--chrom"))
		chromFn = cmdOpts.getOpt("--chrom");

	if(cmdOpts.hasOpt("-s"))
		seqFn = cmdOpts.getOpt("-s");
	if(cmdOpts.hasOpt("--seq"))
		seqFn = cmdOpts.getOpt("--seq");


	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	/* set dbName */
	mtgFn = dbName + METAGENOME_FILE_SUFFIX;
	fmdidxFn = dbName + FMDINDEX_FILE_SUFFIX;

	/* open inputs */
	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	fmdidxIn.open(fmdidxFn.c_str(), ios_base::binary);
	if(!fmdidxIn.is_open()) {
		cerr << "Unable to open '" << fmdidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
	if(!genomeFn.empty()) {
		genomeOut.open(genomeFn.c_str());
		if(!genomeOut.is_open()) {
			cerr << "Unable to write to '" << genomeFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!chromFn.empty()) {
		chromOut.open(chromFn.c_str());
		if(!chromOut.is_open()) {
			cerr << "Unable to write to '" << chromFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!seqFn.empty()) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(seqFn, GZIP_FILE_SUFFIX))
			seqOut.push(boost::iostreams::gzip_compressor());
		else if(StringUtils::endsWith(seqFn, BZIP2_FILE_SUFFIX))
			seqOut.push(boost::iostreams::bzip2_compressor());
		else { }
#endif
		seqOut.push(boost::iostreams::file_sink(seqFn));
	}

	/* load data */
	MetaGenome mtg;
	FMDIndex fmdidx;
	MetaGenomeAnno mta;

	infoLog << "Loading MetaGenome info ..." << endl;
	loadProgInfo(mtgIn);
	if(!mtgIn.bad())
		mtg.load(mtgIn);
	if(mtgIn.bad()) {
		cerr << "Unable to load MetaGenome: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Loading FMD-index ..." << endl;
	loadProgInfo(fmdidxIn);
	if(!fmdidxIn.bad())
		fmdidx.load(fmdidxIn);
	if(fmdidxIn.bad()) {
		cerr << "Unable to load FMD-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Checking MetaGenome indices ..." << endl;
	for(size_t tid = 0; tid < mtg.numChroms(); ++tid) {
		string tname = mtg.getChromName(tid);
		const Genome::Chrom& chr = mtg.getChrom(tid);
		if(!(tname == chr.name && mtg.getChromLoc(tid).length() == chr.size() + 1)) {
			cerr << "Unmatched chrom record for tid: " << tid << endl;
			return EXIT_FAILURE;
		}
	}

	cout << "MetaGenome info: # of genomes: " << mtg.numGenomes() << " size: " << mtg.BDSize() << endl;
	cout << "FMD-index info: length: " << fmdidx.length() << endl;
	cout << "Basic base count:"
			<< " A: " << fmdidx.getBaseCount(DNAalphabet::A)
			<< " C: " << fmdidx.getBaseCount(DNAalphabet::C)
			<< " G: " << fmdidx.getBaseCount(DNAalphabet::G)
			<< " T: " << fmdidx.getBaseCount(DNAalphabet::T)
			<< endl;

	/* if genomeOut is requested */
	if(genomeOut.is_open()) {
		infoLog << "Writing genome list" << endl;
		genomeOut << GENOME_TSV_HEADER << endl;
		for(const Genome& genome : mtg.getGenomes())
			genomeOut << genome.getId() << "\t" << genome.getName() << genome.size() << endl;
	}

	/* if chromOut is requested */
	if(chromOut.is_open()) {
		infoLog << "Writing chrom list" << endl;
		chromOut << CHROM_TSV_HEADER << endl;
		for(const Genome& genome : mtg.getGenomes())
			for(const Genome::Chrom& chr : genome.getChroms())
				chromOut << chr.name << "\t" << chr.size() << "\t" << genome.getId() << "\t" << genome.getName() << endl;
	}

	/* if seqOut requested */
	if(seqOut.is_complete()) {
		infoLog << "Writing chrom sequences" << endl;
		SeqIO seqO(&seqOut, SeqIO::FASTA);
		for(const Genome& genome : mtg.getGenomes())
			for(const Genome::Chrom& chr : genome.getChroms())
				seqO.writeSeq(PrimarySeq(chr.seq, chr.name,
						"genomeId=" + genome.getId() + ";genomeName=" + genome.getName() + ";chromName=" + chr.name));
	}
}
