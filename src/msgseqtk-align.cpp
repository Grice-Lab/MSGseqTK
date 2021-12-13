/*******************************************************************************
 * This file is part of MSGseqTK, a Metagenomics Shot-Gun sequencing ToolKit
 * for ultra-fast and accurate MSG-seq cleaning, mapping and more,
 * based on space-efficient FM-index on entire collection of meta-genomics sequences.
 * Copyright (C) 2020  Qi Zheng
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
 * along with MSGseqTK.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/*
 * msgseqtk-align.cpp
 *
 *  Created on: Dec 11, 2018
 *      Author: zhengqi
 * BAM record exported from this program will have additional application-specified SAM tags:
 * X?: alignment tags, Z? total tags
 * Tag  Type  Description
 * NH   i     Number of reported alignments
 * XN   i     Number of total alignments (primary + secondary) satisfying the user-specified criteria
 * XH   f     alignment log10-likelihood given bases and qualities
 * XP   f     alignment posterior probability, Baysian inferred from XH and prior proportional to XA
 */

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
#include "MSGseqTK.h"
#include "EGUtil.h"
#include "EGSAMtools.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::SAMtools;
using namespace EGriceLab::MSGseqTK;

/* program default values */
//static const int DEFAULT_STRAND = 3;
static const int DEFAULT_NUM_THREADS = 1;
static const int DEFAULT_MAX_REPORT = 1;
static const double DEFAULT_MIN_SCORE_RATE = 0.65;
static const double DEFAULT_MAX_EVALUE = 0.01;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Map/align Metagenomics Shot-Gun (MSG) sequencing reads to database,"
		 << " based on sampling-based Maximal Exact Matched Seeds (MEMS) searches"
		 << " followed by banded Smith-Waterman Dymamic Programming" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif

	cerr << "Usage:    " << progName << "  <DB> <READ-FILE> [MATE-FILE] <-o ALIGN-OUT> [other-options]" << endl
		 << "DB  STR                          : database name/prefix" << endl
		 << "READ-FILE  FILE                  : single-end/forward MSG read file" << ZLIB_SUPPORT << endl
		 << "MATE-FILE  FILE                  : mate/reverse MSG read file" << ZLIB_SUPPORT << endl
		 << "Options:    -o  FILE             : BAM or SAM output file" << endl
		 << "            --global  FLAG       : use Needleman-Wunsch global DP algorithm, no clipping allowed [" << (Alignment::DEFAULT_MODE == Alignment::GLOBAL ? "on" : "off") << "]" << endl
		 << "            --local  FLAG        : use Smith-Waterman local DP algorithm, ends may be soft-clipped [" << (Alignment::DEFAULT_MODE == Alignment::LOCAL ? "on" : "off") << "]" << endl
		 << "            --match  DBL         : score for matches [" << ScoreScheme::DEFAULT_MATCH_SCORE << "]" << endl
		 << "            --mis-match  DBL     : penalty for mis-matches [" << ScoreScheme::DEFAULT_MISMATCH_PENALTY << "]" << endl
		 << "            --gap-open  DBL      : penalty for (affine) gap opening [" << ScoreScheme::DEFAULT_GAP_OPEN_PENALTY << "]" << endl
		 << "            --gap-ext  DBL       : penalty for (affine) gap extension [" << ScoreScheme::DEFAULT_GAP_EXT_PENALTY << "]" << endl
		 << "            -k/--max-report  INT : maximum loci to consider for a read/pair, set to 0 to report all candidate alignments [" << DEFAULT_MAX_REPORT << "]" << endl
		 << "            -s/--min-score       : minimum score rate as a fraction of read-length * match-score [" << DEFAULT_MIN_SCORE_RATE << "]" << endl
		 << "Paired-end:" << endl
		 << "            -m/--mean-ins  DBL   : mean insert size, set to 0 to ignore pairing probabilities (uniform prior) [" << PairingScheme::DEFAULT_MEAN_INSERT << "]" << endl
		 << "            -s/--sd-ins  DBL     : standard deviation of insert size [" << PairingScheme::DEFAULT_CV_INSERT << " * mean]" << endl
		 << "            -I/--min-ins  DBL    : minumum insert size [mean - " << PairingScheme::DEFAULT_OUTLIER_DIST << " * sd]" << endl
		 << "            -X/--max-ins  DBL    : maximum insert size [mean + " << PairingScheme::DEFAULT_OUTLIER_DIST << " * sd]" << endl
		 << "            --no-mixed  FLAG     : suppress unpaired alignments for paired reads" << endl
		 << "            --no-discordant FLAG : suppress discordant alignments for paired reads" << endl
		 << "            --no-tail-over FLAG  : not concordant when mates extend (tail) past each other" << endl
		 << "            --no-contain  FLAG   : not concordant when one mate alignment contains other" << endl
		 << "            --no-overlap  FLAG   : not concordant when mates overlap" << endl
		 << "            --max-npair  INT     : maximum number of pairs to for each forward/reverse mates, 0 for no limit [" << AlignmentPE::MAX_NPAIR << "]" << endl
		 << "Other:" << endl
		 << "            --min-seed  INT      : minimum length of an SMEM to be used as a seed [" << SMEM_LIST::MIN_LENGTH << "]" << endl
		 << "            --max-seed  INT      : maximum length of an SMEM that will trigger re-seeding to avoid missed seeds, 0 for no-reseeding [" << SMEM_LIST::MAX_LENGTH << "]" << endl
		 << "            --max-evalue  DBL    : maximum evalue of an SMEM to be used as a seed, 0 for no limit [" << DEFAULT_MAX_EVALUE << "]" << endl
		 << "            --max-nseed  INT     : maximum # of loci to check for each SMEM [" << SMEM::MAX_NSEED << "]" << endl
//		 << "            --max-lod10  DBL     : maximum log10-liklihood odd (lod10) allowed for a good seed-chain compared to the best seed-chain [" << DEFAULT_MAX_LOD10 << "]" << endl
#ifdef _OPENMP
		 << "            -p|--process INT     : number of threads/cpus for parallel processing [" << DEFAULT_NUM_THREADS << "]" << endl
#endif
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

/**
 * main function to process single-ended reads
 * @return # of reads processed
 */
uint64_t main_SE(const MetaGenome& mtg, const FMDIndex& fmdidx, SeqIO& seqI, SAMfile& out,
		int64_t minSeed, int64_t maxSeed, double maxEvalue, int64_t maxNSeed,
		Alignment::MODE alnMode, double minScoreRate, uint32_t maxReport);

/**
 * main function to process paired-ended reads
 * @return # of read pairs processed
 */
uint64_t main_PE(const MetaGenome& mtg, const FMDIndex& fmdidx, SeqIO& fwdI, SeqIO& revI, SAMfile& out,
		int64_t minSeed, int64_t maxSeed, double maxEvalue, int64_t maxNSeed,
		Alignment::MODE alnMode, double minScoreRate, uint32_t maxReport,
		bool noMixed, bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap, int64_t maxNPair);

/**
 * report and output evaluated Alignments
 * @return # of alignments written
 */
int output(const ALIGN_LIST& alnList, SAMfile& out, uint32_t maxReport, uint16_t extrFlag = 0);

/**
 * report and output evaluated Alignment pairs
 * @return # of pairs written
 */
int output(const PAIR_LIST& pairList, SAMfile& out, uint32_t maxReport);

/** report and output unmapped reads */
int output(const PrimarySeq& read, SAMfile& out);

/** report and output unmapped pairs */
int output(const PrimarySeq& fwdRead, const PrimarySeq& revRead, SAMfile& out);

int main(int argc, char* argv[]) {
	/* variable declarations */
	string fwdInFn, revInFn, outFn;
	string db;
	boost::iostreams::filtering_istream fwdIn, revIn;

	Alignment::MODE alnMode = Alignment::DEFAULT_MODE;
	bool isPaired = false;
	int64_t minSeed = SMEM_LIST::MIN_LENGTH;
	int64_t maxSeed = SMEM_LIST::MAX_LENGTH;
	double maxEvalue = DEFAULT_MAX_EVALUE;
	int64_t maxNSeed = SMEM::MAX_NSEED;
//	double maxLod10 = DEFAULT_MAX_LOD10;
//	double maxIndelRate = DEFAULT_INDEL_RATE;
	int nThreads = DEFAULT_NUM_THREADS;

//	uint32_t maxSeeds = Alignment::MAX_ALIGN;
	uint32_t maxReport = DEFAULT_MAX_REPORT;
	double minScoreRate = DEFAULT_MIN_SCORE_RATE;

	/* pairing options */
	bool noMixed = false;
	bool noDiscordant = false;
	bool noTailOver = false;
	bool noContain = false;
	bool noOverlap = false;
	int64_t maxNPair = AlignmentPE::MAX_NPAIR;

	/* parse options */
	/* main options */
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

	if(!(1 <= cmdOpts.numMainOpts() && cmdOpts.numMainOpts() <= 3)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	db = cmdOpts.getMainOpt(0);
	fwdInFn = cmdOpts.getMainOpt(1);
	if(cmdOpts.numMainOpts() == 3) {
		revInFn = cmdOpts.getMainOpt(2);
		isPaired = true;
	}

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("--global"))
		alnMode = Alignment::GLOBAL;

	if(cmdOpts.hasOpt("--local"))
		alnMode = Alignment::LOCAL;

	if(cmdOpts.hasOpt("--match")) {
		double matchScore = ::atof(cmdOpts.getOptStr("--match"));
		if(!(matchScore >= 0)) {
			cerr << "--match must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		Alignment::ss.setMatchScore(matchScore);
	}

	if(cmdOpts.hasOpt("--mis-match")) {
		double misPenalty = ::atof(cmdOpts.getOptStr("--mis-match"));
		if(!(misPenalty >= 0)) {
			cerr << "--mis-match must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		Alignment::ss.setMismatchPenalty(misPenalty);
	}

	if(cmdOpts.hasOpt("--gap-open")) {
		double op = ::atof(cmdOpts.getOptStr("--gap-open"));
		if(!(op >= 0)) {
			cerr << "--gap-open must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		Alignment::ss.setGapOPenalty(op);
	}

	if(cmdOpts.hasOpt("--gap-ext")) {
		double op = ::atof(cmdOpts.getOptStr("--gap-ext"));
		if(!(op >= 0)) {
			cerr << "--gap-ext must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		Alignment::ss.setGapEPenalty(op);
	}
//
//	if(cmdOpts.hasOpt("--max-seeds"))
//		maxSeeds = ::atoi(cmdOpts.getOptStr("--max-seeds"));

	if(cmdOpts.hasOpt("-k"))
		maxReport = ::atoi(cmdOpts.getOptStr("-k"));
	if(cmdOpts.hasOpt("--max-report"))
		maxReport = ::atoi(cmdOpts.getOptStr("--max-report"));

	if(cmdOpts.hasOpt("-s"))
		minScoreRate = ::atof(cmdOpts.getOptStr("-s"));
	if(cmdOpts.hasOpt("--min-score"))
		minScoreRate = ::atof(cmdOpts.getOptStr("--min-score"));

	/* paired-end options */
	if(cmdOpts.hasOpt({"-m", "--mean-ins"}))
	{
		double meanIns = 0;
		if(cmdOpts.hasOpt("-m"))
			meanIns = ::atof(cmdOpts.getOptStr("-m"));
		if(cmdOpts.hasOpt("--mean-ins"))
			meanIns = ::atof(cmdOpts.getOptStr("--mean-ins"));
		if(!(meanIns >= 0)) {
			cerr << "-m/--mean-ins must non-netative" << endl;
			return EXIT_FAILURE;
		}
		AlignmentPE::ps.setMean(meanIns);
	}

	if(cmdOpts.hasOpt({"-s", "--sd-ins"})) {
		double sdIns = 0;
		if(cmdOpts.hasOpt("-s"))
			sdIns = ::atof(cmdOpts.getOptStr("-s"));
		if(cmdOpts.hasOpt("--sd-ins"))
			sdIns = ::atof(cmdOpts.getOptStr("--sd-ins"));
		if(!(sdIns > 0)) {
			cerr << "-s/--sd-ins must be positive" << endl;
			return EXIT_FAILURE;
		}
		AlignmentPE::ps.setSD(sdIns);
	}
	AlignmentPE::ps.updateRange(); // update range using m and s

	/* set customized ranges */
	if(cmdOpts.hasOpt({"-I", "--min-ins"})) {
		double minIns = 0;
		if(cmdOpts.hasOpt("-I"))
			minIns = ::atof(cmdOpts.getOptStr("-I"));
		if(cmdOpts.hasOpt("--min-ins"))
			minIns = ::atof(cmdOpts.getOptStr("--min-ins"));
		if(!(minIns >= 0)) {
			cerr << "-I/--min-ins must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		AlignmentPE::ps.setMin(minIns);
	}

	if(cmdOpts.hasOpt({"-X", "--max-ins"})) {
		double maxIns = 0;
		if(cmdOpts.hasOpt("-X"))
			maxIns = ::atof(cmdOpts.getOptStr("-X"));
		if(cmdOpts.hasOpt("--max-ins"))
			maxIns = ::atof(cmdOpts.getOptStr("--max-ins"));
		if(!(maxIns > 0)) {
			cerr << "-X/--max-ins must be positive" << endl;
			return EXIT_FAILURE;
		}
		AlignmentPE::ps.setMax(maxIns);
	}

	if(cmdOpts.hasOpt("--no-mixed"))
		noMixed = true;
	if(cmdOpts.hasOpt("--no-discordant"))
		noDiscordant = true;
	if(cmdOpts.hasOpt("--no-tail-over"))
		noTailOver = true;
	if(cmdOpts.hasOpt("--no-contain"))
		noContain = true;
	if(cmdOpts.hasOpt("--no-overlap"))
		noOverlap = true;
	if(cmdOpts.hasOpt("--max-npair"))
		maxNPair = ::atol(cmdOpts.getOptStr("--max-npair"));

	/* other options */
	if(cmdOpts.hasOpt("--min-seed"))
		minSeed = ::atol(cmdOpts.getOptStr("--min-seed"));

	if(cmdOpts.hasOpt("--max-seed"))
		maxSeed = ::atol(cmdOpts.getOptStr("--max-seed"));

	if(cmdOpts.hasOpt("--max-evalue"))
		maxEvalue = ::atof(cmdOpts.getOptStr("--max-evalue"));

	if(cmdOpts.hasOpt("--max-nseed"))
		maxNSeed = ::atol(cmdOpts.getOptStr("--max-nseed"));
//
//	if(cmdOpts.hasOpt("--max-lod10"))
//		maxLod10 = ::atof(cmdOpts.getOptStr("--max-lod10"));

#ifdef _OPENMP
	if(cmdOpts.hasOpt("-p"))
		nThreads = ::atoi(cmdOpts.getOptStr("-p"));
	if(cmdOpts.hasOpt("--process"))
		nThreads = ::atoi(cmdOpts.getOptStr("--process"));
#endif

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(outFn.empty()) {
		cerr << "-o must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(!(maxReport >= 0)) {
		cerr << "--max-report must be non-negagive" << endl;
		return EXIT_FAILURE;
	}

	if(!(minScoreRate <= 1)) {
		cerr << "-s|-min-score must no greater than 1" << endl;
		return EXIT_FAILURE;
	}

	if(!(0 < minSeed)) {
		cerr << "--min-seed must be non-negative" << endl;
		return EXIT_FAILURE;
	}

	if(!(0 == maxSeed || minSeed <= maxSeed)) {
		cerr << "--max-seed must be no smaller than --min-seed" << endl;
		return EXIT_FAILURE;
	}

	if(!(maxEvalue >= 0)) {
		cerr << "--max-evalue must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(maxEvalue == 0)
		maxEvalue = inf;

	if(!(maxNSeed > 0)) {
		cerr << "--max-nseed must be positive" << endl;
		return EXIT_FAILURE;
	}

#ifdef _OPENMP
	if(!(nThreads > 0)) {
		cerr << "-p|--process must be positive" << endl;
		return EXIT_FAILURE;
	}
	omp_set_num_threads(nThreads);
#endif

	/* update set pairing scheme */
	AlignmentPE::ps.updateParam();

	/* guess input seq format */
	SeqIO::FORMAT fmt = SeqIO::guessFormat(fwdInFn);
	if(fmt == SeqIO::UNK) {
		cerr << "Unrecognized sequence format for file '" << fwdInFn << "'" << endl;
		return EXIT_FAILURE;
	}
	bool fixQual = fmt != SeqIO::FASTQ; // non-fastq format need quality fix

	/* open inputs */
#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(fwdInFn, GZIP_FILE_SUFFIX))
		fwdIn.push(boost::iostreams::gzip_decompressor());
	else if(StringUtils::endsWith(fwdInFn, BZIP2_FILE_SUFFIX))
		fwdIn.push(boost::iostreams::bzip2_decompressor());
	else { }
#endif
	boost::iostreams::file_source fwdSrc(fwdInFn);
	if(!fwdSrc.is_open()) {
		cerr << "Unable to open '" << fwdInFn << "' " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	fwdIn.push(fwdSrc);

	if(isPaired) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(revInFn, GZIP_FILE_SUFFIX))
			revIn.push(boost::iostreams::gzip_decompressor());
		else if(StringUtils::endsWith(revInFn, BZIP2_FILE_SUFFIX))
			revIn.push(boost::iostreams::bzip2_decompressor());
		else { }
#endif
		boost::iostreams::file_source revSrc(revInFn);
		if(!revSrc.is_open()) {
			cerr << "Unable to open '" << revInFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		revIn.push(revSrc);
	}

	/* open SeqIO input */
	SeqIO fwdI(dynamic_cast<istream*>(&fwdIn), fmt, fixQual);
	SeqIO revI(dynamic_cast<istream*>(&revIn), fmt, fixQual);
	string mtgFn = db + METAGENOME_FILE_SUFFIX;
	string mgsFn = db + METAGENOME_SEQ_FILE_SUFFIX;
	string fmdidxFn = db + FMDINDEX_FILE_SUFFIX;

	ifstream mtgIn, mgsIn;
	ifstream fmdidxIn;

	MetaGenome mtg;
	FMDIndex fmdidx;

	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	mgsIn.open(mgsFn.c_str(), ios_base::binary);
	if(!mgsIn.is_open()) {
		cerr << "Unable to open '" << mgsFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	fmdidxIn.open(fmdidxFn.c_str(), ios_base::binary);
	if(!fmdidxIn.is_open()) {
		cerr << "Unable to open '" << fmdidxFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
	string mode;
	if(StringUtils::endsWith(outFn, SAM_SUFFIX))
		mode = "w";
	else if(StringUtils::endsWith(outFn, BAM_SUFFIX))
		mode = "wb";
	else {
		cerr << "Unable to determine SAM/BAM output mode" << endl;
		return EXIT_FAILURE;
	}
	SAMfile out(outFn, mode);

	/* load data */
	infoLog << "Loading MetaGenome info ..." << endl;
	loadProgInfo(mtgIn);
	if(!mtgIn.bad())
		mtg.load(mtgIn);
	if(mtgIn.bad()) {
		cerr << "Unable to load MetaGenome info: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Loading MetaGenome seq ..." << endl;
	mtg.loadSeq(mgsIn);
	if(mgsIn.bad()) {
		cerr << "Unable to load MetaGenome seq: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Loading FMD-index ..." << endl;
	loadProgInfo(fmdidxIn);
	if(!fmdidxIn.bad())
		fmdidx.load(fmdidxIn);
	if(fmdidxIn.bad()) {
		cerr << "Unable to load reference FMD-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	map<string, uint32_t> targetLen;
	for(size_t i = 0; i < mtg.numChroms(); ++i)
		targetLen[mtg.getChromName(i)] = mtg.getChrom(i).size();

	BAMheader header(mtg.getChromNames(), targetLen);

	/* add header and additional tags */
	header.setHDTag("SQ", "unsorted");
	header.addTag("@PG", "ID:" + progName + " VN:" + progVer.toString() + " PN:" + cmdOpts.getProg() + " CL:" + cmdOpts.getCmdStr());
	out.setHeader(header);
	out.writeHeader();

	/* main processing */
	infoLog << "Aligning input reads ..." << endl;
	uint64_t nProcessed = 0;
	chrono::time_point<chrono::steady_clock> start = chrono::steady_clock::now();
	if(!isPaired)
		nProcessed = main_SE(mtg, fmdidx, fwdI, out, minSeed, maxSeed, maxEvalue, maxNSeed, alnMode, minScoreRate, maxReport);
	else
		nProcessed = main_PE(mtg, fmdidx, fwdI, revI, out, minSeed, maxSeed, maxEvalue, maxNSeed, alnMode, minScoreRate, maxReport,
				noMixed, noDiscordant, noTailOver, noContain, noOverlap, maxNPair);
	chrono::time_point<chrono::steady_clock> fin = chrono::steady_clock::now();
	infoLog << "Read alignment finished. Total processed reads: " <<  nProcessed
			<< " . Elapsed time: "
			<< chrono::duration_cast<std::chrono::seconds>(fin - start).count()
			<< " sec" << endl;
}

int output(const ALIGN_LIST& alnList, SAMfile& out, uint32_t maxReport, uint16_t extraFlag) {
	/* export BAM records */
	size_t numReport = maxReport == 0 ? alnList.size() : std::min((size_t) maxReport, alnList.size());
	/* generate BAM records first, which does not need locks or openMP critical pragma */
	for(size_t i = 0; i < numReport; ++i) {
		const Alignment& aln = alnList[i];
		/* add a newly constructed BAM */
		BAM bamAln = aln.exportBAM();
		/* set additional flags only available at output */
		bamAln.setSecondaryFlag(i > 0);
		bamAln.setFlag(bamAln.getFlag() | extraFlag);
		/* set standard aux tags */
		bamAln.setAux(Alignment::NUM_REPORTED_ALIGNMENT_TAG, numReport);
		bamAln.setAux(Alignment::NUM_TOTAL_ALIGNMENT_TAG, alnList.size());
		/* set customized aux tags */
		bamAln.setAux(Alignment::ALIGNMENT_POSTERIOR_PROB_TAG, aln.getPostP());
		out.write(bamAln);
	}
	return numReport;
}

int output(const PAIR_LIST& pairList, SAMfile& out, uint32_t maxReport) {
	/* export BAM records */
	size_t numReport = maxReport == 0 ? pairList.size() : std::min<size_t>(maxReport, pairList.size());
	/* generate BAM records lists first, which does not need locks or openMP critical */
	for(size_t i = 0; i < numReport; ++i) {
		const AlignmentPE& pair = pairList[i];
		/* add newly constructed BAM */
		BAM fwdBam = pair.exportFwdBAM();
		BAM revBam = pair.exportRevBAM();
		/* set additional flags only available at output */
		fwdBam.setSecondaryFlag(i > 0);
		revBam.setSecondaryFlag(i > 0);
		/* set standard aux tags */
		fwdBam.setAux(Alignment::NUM_REPORTED_ALIGNMENT_TAG, numReport);
		fwdBam.setAux(Alignment::NUM_TOTAL_ALIGNMENT_TAG, pairList.size());
		revBam.setAux(Alignment::NUM_REPORTED_ALIGNMENT_TAG, numReport);
		revBam.setAux(Alignment::NUM_TOTAL_ALIGNMENT_TAG, pairList.size());
		/* set customized aux tags */
		fwdBam.setAux(Alignment::ALIGNMENT_POSTERIOR_PROB_TAG, pair.postP);
		revBam.setAux(Alignment::ALIGNMENT_POSTERIOR_PROB_TAG, pair.postP);
		out.write(fwdBam);
		out.write(revBam);
	}
	return numReport;
}

int output(const PrimarySeq& read, SAMfile& out) {
	return out.write(BAM(read.getName(), read.length(), dna::nt16Encode(read.getSeq()), read.getQual(), BAM_FUNMAP));
}

int output(const PrimarySeq& fwdRead, const PrimarySeq& revRead, SAMfile& out) {
	out.write(BAM(fwdRead.getName(), fwdRead.length(), dna::nt16Encode(fwdRead.getSeq()), fwdRead.getQual(), BAM_FUNMAP | BAM_FMUNMAP | BAM_FPAIRED | BAM_FREAD1));
	out.write(BAM(revRead.getName(), revRead.length(), dna::nt16Encode(revRead.getSeq()), revRead.getQual(), BAM_FUNMAP | BAM_FMUNMAP | BAM_FPAIRED | BAM_FREAD2));
	return 2;
}

uint64_t main_SE(const MetaGenome& mtg, const FMDIndex& fmdidx, SeqIO& seqI, SAMfile& out,
		int64_t minSeed, int64_t maxSeed, double maxEvalue, int64_t maxNSeed,
		Alignment::MODE alnMode, double minScoreRate, uint32_t maxReport) {
	uint64_t nRead = 0;
	/* search MEMS for each read */
#pragma omp parallel
	{
#pragma omp master
		{
			while(seqI.hasNext()) {
				const PrimarySeq read = seqI.nextSeq();
				const PrimarySeq rcRead = read.revcom();
#pragma omp task firstprivate(read, rcRead)
				{
					/* get SeedList */
					SeedList seeds = SMEM_LIST::findSeeds(&read, &mtg, &fmdidx, minSeed, maxSeed, maxEvalue, maxNSeed);
					if(seeds.empty()) {
#pragma omp critical(LOG)
						debugLog << "Unable to find any valid SMEM seads for '" << read.getName() << "', ignore" << endl;
#pragma omp critical(BAM_OUTPUT)
						output(read, out);
					}
					else {
						const int64_t L = read.length();
						const int64_t maxMismatch = std::ceil(L * Alignment::MAX_MISMATCH_RATE);
						const int64_t maxIndel = std::ceil(L * Alignment::MAX_INDEL_RATE);
						/* get SeedChains */
						ChainList chains = SeedChain::getChains(seeds, maxMismatch, maxIndel);
						/* filter chains by log10-odd (lod10) */
//						SeedChain::filter(chains, maxLod10);
						/* get unique chains */
						SeedChain::uniq(chains);
						/* get alignments from SeedMatchList */
						ALIGN_LIST alnList = Alignment::buildAlignments(&read, &rcRead, mtg, chains, alnMode);
						/* filter alignments */
						Alignment::filter(alnList, minScoreRate);
						if(alnList.empty()) { // no good alignment found
#pragma omp critical(BAM_OUTPUT)
							output(read, out);
						}
						else {
							/* evaluate alignments */
							Alignment::evaluate(alnList);
							/* calculate mapQ for alignments */
							Alignment::calcMapQ(alnList);
							/* sort alignments */
							Alignment::sort(alnList);
							/* output alignments */
#pragma omp critical(BAM_OUTPUT)
							output(alnList, out, maxReport);
						}
					} /* end task */
#pragma omp atomic
					nRead++;
				} /* end each read */
			} /* end master, implicit barrier */
		} /* end parallel */
	}
	return nRead;
}

uint64_t main_PE(const MetaGenome& mtg, const FMDIndex& fmdidx, SeqIO& fwdI, SeqIO& revI, SAMfile& out,
		int64_t minSeed, int64_t maxSeed, double maxEvalue, int64_t maxNSeed,
		Alignment::MODE alnMode, double minScoreRate, uint32_t maxReport,
		bool noMixed, bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap, int64_t maxNPair) {
	uint64_t nPair = 0;
	/* search MEMS for each read */
#pragma omp parallel
	{
#pragma omp master 
		{
			while(fwdI.hasNext() && revI.hasNext()) {
				PrimarySeq fwdRead = fwdI.nextSeq();
				PrimarySeq revRead = revI.nextSeq();

				if(fwdRead.getName() != revRead.getName()) {
					fwdRead.trimNameExt();
					revRead.trimNameExt();
					if(fwdRead.getName() != revRead.getName()) { // error
						cerr << "Error: unmatched read ID between " << fwdRead.getName() << " and " << revRead.getName() << endl;
						abort();
					}
				}
				PrimarySeq rcFwdRead = static_cast<const PrimarySeq&>(fwdRead).revcom();
				PrimarySeq rcRevRead = static_cast<const PrimarySeq&>(revRead).revcom();
#pragma omp task firstprivate(fwdRead, revRead, rcFwdRead, rcRevRead)
				{
					SeedListPE seedsPE = SMEM_LIST::findSeedsPE(&fwdRead, &revRead, &mtg, &fmdidx, minSeed, maxSeed, maxEvalue, maxNSeed);
					if(seedsPE.first.empty() && seedsPE.second.empty()) {
#pragma omp critical(LOG)
						debugLog << "Unable to find any valid SMEMS seeds for read pair '" << fwdRead.getName() << "'" << endl;
#pragma omp critical(BAM_OUTPUT)
						output(fwdRead, revRead, out);
					}
					else {
						const int64_t fwdL = fwdRead.length();
						const int64_t revL = revRead.length();
						const int64_t fwdMaxMismatch = std::ceil(fwdL * Alignment::MAX_MISMATCH_RATE);
						const int64_t fwdMaxIndel = std::ceil(fwdL * Alignment::MAX_INDEL_RATE);
						const int64_t revMaxMismatch = std::ceil(revL * Alignment::MAX_MISMATCH_RATE);
						const int64_t revMaxIndel = std::ceil(revL * Alignment::MAX_INDEL_RATE);
						/* get chains from seeds */
						ChainList fwdChains = SeedChain::getChains(seedsPE.first, fwdMaxMismatch, fwdMaxIndel);
						ChainList revChains = SeedChain::getChains(seedsPE.second, revMaxMismatch, revMaxIndel);
						/* filter chains */
//						SeedChain::filter(fwdChains, maxLod10);
//						SeedChain::filter(revChains, maxLod10);
						/* get unique chains */
						SeedChain::uniq(fwdChains);
						SeedChain::uniq(revChains);
						/* get alignments */
						ALIGN_LIST fwdAlnList = Alignment::buildAlignments(&fwdRead, &rcFwdRead, mtg, fwdChains, alnMode);
						ALIGN_LIST revAlnList = Alignment::buildAlignments(&revRead, &rcRevRead, mtg, revChains, alnMode);
						/* 1st pass filtes by alignment scores */
						Alignment::filter(fwdAlnList, minScoreRate);
						Alignment::filter(revAlnList, minScoreRate);
						if(fwdAlnList.empty() && revAlnList.empty()) { // both mate cannot map
#pragma omp critical(BAM_OUTPUT)
							output(fwdRead, revRead, out);
						}
						else {
							/* evaluate alignments */
							Alignment::evaluate(fwdAlnList);
							Alignment::evaluate(revAlnList);
							/* get pairs */
							PAIR_LIST pairList = AlignmentPE::getPairs(fwdAlnList, revAlnList);
							/* sort pairs by loglik decreasingly */
							AlignmentPE::sort(pairList);
							/* filter pairs */
							AlignmentPE::filter(pairList, noDiscordant, noTailOver, noContain, noOverlap, maxNPair);
							if(!pairList.empty()) {
								/* calculate mapQ for pairs */
								AlignmentPE::calcMapQ(pairList);
								/* output pairs */
#pragma omp critical(BAM_OUTPUT)
								output(pairList, out, maxReport);
							}
							else if(!noMixed) { /* pair-end matching failed, try unpaired */
								if(!fwdAlnList.empty()) {
									//								debugLog << "Alignment pairing failed, reporting unpaired alignment for forward read " << fwdRead.getName() << endl;
									/* calculate mapQ */
									Alignment::calcMapQ(fwdAlnList);
									/* sort alignments */
									Alignment::sort(fwdAlnList);
									/* output alignments */
#pragma omp critical(BAM_OUTPUT)
									output(fwdAlnList, out, maxReport, BAM_FPAIRED | BAM_FREAD1);
								}
								if(!revAlnList.empty()) {
									//								debugLog << "Alignment pairing failed, reporting unpaired alignment for reverse read " << revRead.getName() << endl;
									/* calculate mapQ */
									Alignment::calcMapQ(revAlnList);
									/* sort alignments */
									Alignment::sort(revAlnList);
									/* output alignments */
#pragma omp critical(BAM_OUTPUT)
									output(revAlnList, out, maxReport, BAM_FPAIRED | BAM_FREAD2);
								}
							}
						}
					} /* end SeedMatch tests */
				} /* end task */
#pragma omp atomic
				nPair++;
			} /* end each read */
		} /* end master, implicit barrier */
	} /* end parallel */
	return nPair;
}
