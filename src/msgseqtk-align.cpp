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
static const double DEFAULT_BEST_FRAC = 0.85;
static const int DEFAULT_MIN_INSERT = 0;
static const int DEFAULT_MAX_INSERT = 750;
static const string NUM_REPORTED_ALIGNMENT_TAG = "NH";
static const string MISMATCH_POSITION_TAG = "MD";
static const string NUM_TOTAL_ALIGNMENT_TAG = "XN";
static const string ALIGNMENT_LOG10LIK_TAG = "XH";
static const string ALIGNMENT_POSTERIOR_PROB_TAG = "XP";

/* static object for this program */
static ScoreScheme ss;

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
		 << "            --match  DBL         : score for matches [" << ScoreScheme::DEFAULT_MATCH_SCORE << "]" << endl
		 << "            --mis-match  DBL     : penalty for mis-matches [" << ScoreScheme::DEFAULT_MISMATCH_PENALTY << "]" << endl
		 << "            --gap-open  DBL      : penalty for (affine) gap opening [" << ScoreScheme::DEFAULT_GAP_OPEN_PENALTY << "]" << endl
		 << "            --gap-ext  DBL       : penalty for (affine) gap extension [" << ScoreScheme::DEFAULT_GAP_EXT_PENALTY << "]" << endl
		 << "            --max-mems  INT      : maximum # of different loc/MEMS to check for a read/pair [" << Alignment::MAX_ALIGN << "]" << endl
		 << "            -k/--max-report  INT : maximum loci to consider for a read/pair, set to 0 to report all candidate alignments [" << DEFAULT_MAX_REPORT << "]" << endl
		 << "            -f--best-frac        : minimum score as a fraction of the highest alignment score of all candidates to consider for full evaluation [" << DEFAULT_BEST_FRAC << "]" << endl
		 << "Paired-end:" << endl
		 << "            -I/--min-ins  INT    : minumum insert size [" << DEFAULT_MIN_INSERT << "]" << endl
		 << "            -X/--max-ins  INT    : maximum insert size, 0 for no limit [" << DEFAULT_MAX_INSERT << "]" << endl
		 << "            --no-mixed  FLAG     : suppress unpaired alignments for paired reads" << endl
		 << "            --no-discordant FLAG : suppress discordant alignments for paired reads" << endl
		 << "            --no-tail-over FLAG  : not concordant when mates extend (tail) past each other" << endl
		 << "            --no-contain  FLAG   : not concordant when one mate alignment contains other" << endl
		 << "            --no-overlap  FLAG   : not concordant when mates overlap" << endl
		 << "Other:" << endl
		 << "            -e|--evalue  DBL     : maximum e-value to consider an MEM as significant [" << MEMS::DEFAULT_MAX_EVALUE << "]" << endl
#ifdef _OPENMP
		 << "            -p|--process INT     : number of threads/cpus for parallel processing [" << DEFAULT_NUM_THREADS << "]" << endl
#endif
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

/**
 * main function to process single-ended reads
 * @return 0 if success, return non-zero otherwise
 */
int main_SE(const MetaGenome& mtg, const FMDIndex& fmdidx,
		SeqIO& seqI, SAMfile& out, double maxEvalue,
		uint32_t maxSeedMatches, double bestFrac, uint32_t maxReport);

/**
 * main function to process paired-ended reads
 * @return 0 if success, return non-zero otherwise
 */
int main_PE(const MetaGenome& mtg, const FMDIndex& fmdidx,
		SeqIO& fwdI, SeqIO& revI, SAMfile& out, double maxEvalue,
		uint32_t maxSeedMatches, double bestFrac, uint32_t maxReport,
		int32_t minIns, int32_t maxIns,
		bool noMixed, bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap);

/**
 * report and output evaluated Alignments
 * @return # of alignments written
 */
int output(const ALIGN_LIST& alnList, SAMfile& out, uint32_t maxReport);

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

	bool isPaired = false;
	double maxEvalue = MEMS::DEFAULT_MAX_EVALUE;
//	double maxIndelRate = DEFAULT_INDEL_RATE;
	int nThreads = DEFAULT_NUM_THREADS;

	uint32_t maxSeedMatches = Alignment::MAX_ALIGN;
	uint32_t maxReport = DEFAULT_MAX_REPORT;
	double bestFrac = DEFAULT_BEST_FRAC;

	/* mate options */
	int32_t minIns = DEFAULT_MIN_INSERT;
	int32_t maxIns = DEFAULT_MAX_INSERT;
	bool noMixed = false;
	bool noDiscordant = false;
	bool noTailOver = false;
	bool noContain = false;
	bool noOverlap = false;

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

	if(cmdOpts.hasOpt("--match")) {
		double matchScore = ::atof(cmdOpts.getOptStr("--match"));
		if(!(matchScore >= 0)) {
			cerr << "--match must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		ss.setMatchScore(matchScore);
	}

	if(cmdOpts.hasOpt("--mis-match")) {
		double misPenalty = ::atof(cmdOpts.getOptStr("--mis-match"));
		if(!(misPenalty >= 0)) {
			cerr << "--mis-match must be non-negative" << endl;
			return EXIT_FAILURE;
		}
		ss.setMismatchPenalty(misPenalty);
	}

	if(cmdOpts.hasOpt("--gap-open")) {
		ss.gapOPenalty = ::atof(cmdOpts.getOptStr("--gap-open"));
		if(!(ss.gapOPenalty >= 0)) {
			cerr << "--gap-open must be non-negative" << endl;
			return EXIT_FAILURE;
		}
	}

	if(cmdOpts.hasOpt("--gap-ext")) {
		ss.gapEPenalty = ::atof(cmdOpts.getOptStr("--gap-ext"));
		if(!(ss.gapEPenalty >= 0)) {
			cerr << "--gap-ext must be non-negative" << endl;
			return EXIT_FAILURE;
		}
	}

	if(cmdOpts.hasOpt("--max-mems"))
		maxSeedMatches = ::atoi(cmdOpts.getOptStr("--max-mems"));

	if(cmdOpts.hasOpt("-k"))
		maxReport = ::atoi(cmdOpts.getOptStr("-k"));
	if(cmdOpts.hasOpt("--max-report"))
		maxReport = ::atoi(cmdOpts.getOptStr("--max-report"));

	if(cmdOpts.hasOpt("--best-frac"))
		bestFrac = ::atof(cmdOpts.getOptStr("--best-frac"));

	/* paired-end options */
	if(cmdOpts.hasOpt("-I"))
		minIns = ::atoi(cmdOpts.getOptStr("-I"));
	if(cmdOpts.hasOpt("--min-ins"))
		minIns = ::atoi(cmdOpts.getOptStr("--min-ins"));

	if(cmdOpts.hasOpt("-X"))
		maxIns = ::atoi(cmdOpts.getOptStr("-X"));
	if(cmdOpts.hasOpt("--max-ins"))
		maxIns = ::atoi(cmdOpts.getOptStr("--max-ins"));

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

	/* other options */
	if(cmdOpts.hasOpt("-e"))
		maxEvalue = ::atof(cmdOpts.getOptStr("-e"));
	if(cmdOpts.hasOpt("--evalue"))
		maxEvalue = ::atoi(cmdOpts.getOptStr("--evalue"));

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

	if(!(maxSeedMatches > 0)) {
		cerr << "--max-mems must be positive" << endl;
		return EXIT_FAILURE;
	}

	if(!(maxReport >= 0)) {
		cerr << "--max-report must be non-negagive" << endl;
		return EXIT_FAILURE;
	}

	if(!(0 <= bestFrac && bestFrac <= 1)) {
		cerr << "--best-frac must between 0 and 1" << endl;
		return EXIT_FAILURE;
	}

	if(!(minIns >= 0)) {
		cerr << "-I/--min-ins must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(!(maxIns >= minIns)) {
		cerr << "-X/--max-ins must be not smaller than -I/--min-ins" << endl;
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
	string fmdidxFn = db + FMDINDEX_FILE_SUFFIX;

	ifstream mtgIn;
	ifstream fmdidxIn;

	MetaGenome mtg;
	FMDIndex fmdidx;

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
		cerr << "Unable to load MetaGenome: " << ::strerror(errno) << endl;
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
	for(size_t i = 0; i < mtg.getChromNames().size(); ++i)
		targetLen[mtg.getChromName(i)] = mtg.getChromLen(i);

	BAMheader header(mtg.getChromNames(), targetLen);

	/* add header and additional tags */
	header.setHDTag("SQ", "unsorted");
	header.addTag("@PG", "ID:" + progName + " VN:" + progVer.toString() + " PN:" + cmdOpts.getProg() + " CL:" + cmdOpts.getCmdStr());
	out.setHeader(header);
	out.writeHeader();

	/* main processing */
	if(!isPaired)
		return main_SE(mtg, fmdidx, fwdI, out, maxEvalue, maxSeedMatches, bestFrac, maxReport);
	else
		return main_PE(mtg, fmdidx, fwdI, revI, out, maxEvalue, maxSeedMatches, bestFrac, maxReport,
				minIns, maxIns, noMixed, noDiscordant, noTailOver, noContain, noOverlap);
}

int output(const ALIGN_LIST& alnList, SAMfile& out, uint32_t maxReport) {
	/* export BAM records */
	size_t numReport = maxReport == 0 ? alnList.size() : std::min((size_t) maxReport, alnList.size());
	for(size_t i = 0; i < numReport; ++i) {
		const Alignment& aln = alnList[i];
		/* construct BAM */
		BAM bamAln = aln.exportBAM();
		/* set flags */
		bamAln.setSecondaryFlag(i > 0);
		/* set standard aux tags */
		bamAln.setAux(NUM_REPORTED_ALIGNMENT_TAG, numReport);
		bamAln.setAux(NUM_TOTAL_ALIGNMENT_TAG, alnList.size());
		bamAln.setAux(MISMATCH_POSITION_TAG, aln.getAlnMDTag());
		/* set customized aux tags */
		bamAln.setAux(ALIGNMENT_LOG10LIK_TAG, aln.log10P);
		bamAln.setAux(ALIGNMENT_POSTERIOR_PROB_TAG, aln.postP);
		/* write BAM */
		out.write(bamAln);
	}
	return numReport;
}

int output(const PAIR_LIST& pairList, SAMfile& out, uint32_t maxReport) {
	/* export BAM records */
	size_t numReport = maxReport == 0 ? pairList.size() : std::min((size_t) maxReport, pairList.size());
	for(size_t i = 0; i < numReport; ++i) {
		const AlignmentPE& pair = pairList[i];
		/* construct BAM */
		BAM bamFwd = pair.exportFwdBAM();
		BAM bamRev = pair.exportRevBAM();
		/* set flags */
		bamFwd.setSecondaryFlag(i > 0);
		bamRev.setSecondaryFlag(i > 0);
		/* set mapQ */
		bamFwd.setMapQ(pair.mapQ);
		bamRev.setMapQ(pair.mapQ);
		/* set standard aux tags */
		bamFwd.setAux(NUM_REPORTED_ALIGNMENT_TAG, numReport);
		bamFwd.setAux(NUM_TOTAL_ALIGNMENT_TAG, pairList.size());
		bamFwd.setAux(MISMATCH_POSITION_TAG, pair.fwdAln->getAlnMDTag());

		bamRev.setAux(NUM_REPORTED_ALIGNMENT_TAG, numReport);
		bamRev.setAux(NUM_TOTAL_ALIGNMENT_TAG, pairList.size());
		bamRev.setAux(MISMATCH_POSITION_TAG, pair.revAln->getAlnMDTag());
		/* set customized aux tags */
		bamFwd.setAux(ALIGNMENT_LOG10LIK_TAG, pair.log10lik());
		bamFwd.setAux(ALIGNMENT_POSTERIOR_PROB_TAG, pair.postP);
		bamRev.setAux(ALIGNMENT_LOG10LIK_TAG, pair.log10lik());
		bamRev.setAux(ALIGNMENT_POSTERIOR_PROB_TAG, pair.postP);
		/* write BAM */
		out.write(bamFwd);
		out.write(bamRev);
	}
	return numReport;
}

int output(const PrimarySeq& read, SAMfile& out) {
	return out.write(BAM(read.getName(), read.length(), dna::nt16Encode(read.getSeq()), read.getQual()));
}

int output(const PrimarySeq& fwdRead, const PrimarySeq& revRead, SAMfile& out) {
	out.write(BAM(fwdRead.getName(), fwdRead.length(), dna::nt16Encode(fwdRead.getSeq()), fwdRead.getQual(), BAM_FREAD1));
	out.write(BAM(revRead.getName(), revRead.length(), dna::nt16Encode(revRead.getSeq()), revRead.getQual(), BAM_FREAD2));
	return 2;
}

int main_SE(const MetaGenome& mtg, const FMDIndex& fmdidx,
		SeqIO& seqI, SAMfile& out, double maxEvalue,
		uint32_t maxSeedMatches, double bestFrac, uint32_t maxReport) {
	infoLog << "Aligning input reads" << endl;
	/* search MEMS for each read */
#pragma omp parallel
	{
#pragma omp master
		{
			while(seqI.hasNext()) {
				const PrimarySeq& read = seqI.nextSeq();
#pragma omp task firstprivate(read)
				{
					MEMS mems = MEMS::searchMEMS(&read, &mtg, &fmdidx, maxEvalue, GLoc::FWD); // fwd sampling
					mems.findLocs(); // only find fwd mapping locs (which can be on reverse-complement genome)
					/* get SeedMatchList */
					const Alignment::SeedMatchList& seedMatches = Alignment::getSeedMatchList(mems, maxSeedMatches);
					if(seedMatches.empty()) {
#pragma omp critical(LOG)
						debugLog << "Unable to find any valid SeedMatches for '" << read.getName() << "'" << endl;
#pragma omp critical(BAM_OUTPUT)
						output(read, out);
					}
					else {
						const PrimarySeq& rcRead = read.revcom();
						debugLog << "id: " << read.getName() << endl;
						/* get alignments from SeedMatchList */
						debugLog << "found " << seedMatches.size() << " seed-matches" << endl;
						ALIGN_LIST alnList = Alignment::getAlignments(&ss, &mtg, &read, &rcRead, seedMatches);
						debugLog << "get " << alnList.size() << " potential aligns" << endl;
						/* filter alignments */
						Alignment::filter(alnList, bestFrac);
						debugLog << "get " << alnList.size() << " best aligns" << endl;
						/* evaluate alignments */
						Alignment::evaluate(alnList);
						for(const Alignment& aln : alnList) {
							debugLog << "cigar:" << BAM::decodeCigar(aln.getAlnCigar())
							<< " aln.loglik10():" << aln.log10lik() << endl;
						}
						/* calculate mapQ for alignments */
						Alignment::calcMapQ(alnList);
						/* sort alignments */
						Alignment::sort(alnList);
						/* output alignments */
#pragma omp critical(BAM_OUTPUT)
						output(alnList, out, maxReport);
					} /* end task */
				} /* end each read */
			} /* end master, implicit barrier */
		} /* end parallel */
	}
	return out.close();
}

int main_PE(const MetaGenome& mtg, const FMDIndex& fmdidx,
		SeqIO& fwdI, SeqIO& revI, SAMfile& out, double maxEvalue,
		uint32_t maxSeedMatches, double bestFrac, uint32_t maxReport,
		int32_t minIns, int32_t maxIns,
		bool noMixed, bool noDiscordant, bool noTailOver, bool noContain, bool noOverlap) {
	infoLog << "Aligning paired-end reads" << endl;
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
#pragma omp task firstprivate(fwdRead, revRead)
				{
					MEMS_PE memsPE = MEMS_PE::sampleMEMS(&fwdRead, &revRead, &mtg, &fmdidx);
					memsPE.findLocs(); // find locs on fwd strand only
					/* get SeedMatchLists */
					Alignment::SeedMatchList fwdSeedMatches = Alignment::getSeedMatchList(memsPE.fwdMems, maxSeedMatches);
					Alignment::SeedMatchList revSeedMatches = Alignment::getSeedMatchList(memsPE.revMems, maxSeedMatches);
					if(fwdSeedMatches.empty() && revSeedMatches.empty()) {
#pragma omp critical(LOG)
						debugLog << "Unable to find any valid SeedMatches for read pair '" << fwdRead.getName() << "'" << endl;
#pragma omp critical(BAM_OUTPUT)
						output(fwdRead, revRead, out);
					}
					else {
						PrimarySeq rcFwdRead = static_cast<const PrimarySeq&>(fwdRead).revcom();
						PrimarySeq rcRevRead = static_cast<const PrimarySeq&>(revRead).revcom();
						ALIGN_LIST fwdAlnList = Alignment::getAlignments(&ss, &mtg, &fwdRead, &rcFwdRead, fwdSeedMatches);
						ALIGN_LIST revAlnList = Alignment::getAlignments(&ss, &mtg, &revRead, &rcRevRead, revSeedMatches);
						/* filter alignments */
						Alignment::filter(fwdAlnList, bestFrac);
						Alignment::filter(revAlnList, bestFrac);
						/* evaluate alignments */
						Alignment::evaluate(fwdAlnList);
						Alignment::evaluate(revAlnList);
						/* get pairs */
						PAIR_LIST pairList = AlignmentPE::getPairs(fwdAlnList, revAlnList);
						/* filter pairs */
						AlignmentPE::filter(pairList, minIns, maxIns, noDiscordant, noTailOver, noContain, noOverlap);
						if(!pairList.empty()) {
							/* calculate mapQ for pairs */
							AlignmentPE::calcMapQ(pairList);
							/* sort pairs */
							AlignmentPE::sort(pairList);
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
								output(fwdAlnList, out, maxReport);
							}
							if(!revAlnList.empty()) {
//								debugLog << "Alignment pairing failed, reporting unpaired alignment for reverse read " << revRead.getName() << endl;
								/* calculate mapQ */
								Alignment::calcMapQ(revAlnList);
								/* sort alignments */
								Alignment::sort(revAlnList);
								/* output alignments */
#pragma omp critical(BAM_OUTPUT)
								output(revAlnList, out, maxReport);
							}
						}
					} /* end SeedMatch tests */
				} /* end task */
			} /* end each read */
		} /* end master, implicit barrier */
	} /* end parallel */
	return out.close();
}
