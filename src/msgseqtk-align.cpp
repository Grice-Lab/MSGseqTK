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
static const int DEFAULT_STRAND = 3;
static const int DEFAULT_NUM_THREADS = 1;
static const int DEFAULT_MAX_REPORT = 1;
static const double DEFAULT_BEST_FRAC = 0.85;
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
		 << "DB  STR                          : database name" << endl
		 << "READ-FILE  FILE                  : single-end/forward MSG read file" << ZLIB_SUPPORT << endl
		 << "MATE-FILE  FILE                  : mate/reverse MSG read file" << ZLIB_SUPPORT << endl
		 << "Options:    -o  FILE             : BAM or SAM output file" << endl
		 << "            --match  DBL         : score for matches [" << ScoreScheme::DEFAULT_MATCH_SCORE << "]" << endl
		 << "            --mis-match  DBL     : penalty for mis-matches [" << ScoreScheme::DEFAULT_MISMATCH_PENALTY << "]" << endl
		 << "            --gap-open  DBL      : penalty for (affine) gap opening [" << ScoreScheme::DEFAULT_GAP_OPEN_PENALTY << "]" << endl
		 << "            --gap-ext  DBL       : penalty for (affine) gap extension [" << ScoreScheme::DEFAULT_GAP_EXT_PENALTY << "]" << endl
		 << "            --max-mems  INT      : maximum # of different loc/MEMS to check for a read/pair [" << Alignment::MAX_ITER << "]" << endl
		 << "            --max-report  INT    : maximum loci to consider for a read/pair, set to 0 to report all candidate alignments [" << DEFAULT_MAX_REPORT << "]" << endl
		 << "            -f--best-frac        : minimum score as a fraction of the best alignment to consider as a candidate for full evaluation [" << DEFAULT_BEST_FRAC << "]" << endl
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
int main_SE(const MetaGenome& mtg, const FMIndex& fmidx,
		SeqIO& seqI, SAMfile& out, RNG& rng, int strand,
		uint32_t maxMems, double bestFrac, uint32_t maxReport);

/**
 * main function to process paired-ended reads
 * @return 0 if success, return non-zero otherwise
 */
int main_PE(const MetaGenome& mtg, const FMIndex& fmidx,
		SeqIO& fwdI, SeqIO& revI, SAMfile& out, RNG& rng, int strand,
		uint32_t maxMems, double bestFrac, uint32_t maxReport);

int main(int argc, char* argv[]) {
	/* variable declarations */
	string fwdInFn, revInFn, outFn;
	string db;
	boost::iostreams::filtering_istream fwdIn, revIn;

	int strand = DEFAULT_STRAND;
	unsigned seed = time(NULL); // using time as default seed
//	double maxIndelRate = DEFAULT_INDEL_RATE;
	int nThreads = DEFAULT_NUM_THREADS;

	typedef boost::random::mt11213b RNG; /* random number generator type */

	uint32_t maxMems = Alignment::MAX_ITER;
	uint32_t maxReport = DEFAULT_MAX_REPORT;
	double bestFrac = DEFAULT_BEST_FRAC;

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

	if(!(1 <= cmdOpts.numMainOpts() && cmdOpts.numMainOpts() <= 3)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	db = cmdOpts.getMainOpt(0);
	fwdInFn = cmdOpts.getMainOpt(1);
	if(cmdOpts.numMainOpts() == 3)
		revInFn = cmdOpts.getMainOpt(2);

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
		maxMems = ::atoi(cmdOpts.getOptStr("--max-mems"));

	if(cmdOpts.hasOpt("--max-report"))
		maxReport = ::atoi(cmdOpts.getOptStr("--max-report"));

	if(cmdOpts.hasOpt("--best-frac"))
		bestFrac = ::atof(cmdOpts.getOptStr("--best-frac"));

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
	if(outFn.empty()) {
		cerr << "-o must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(!(maxMems > 0)) {
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
	string mtgFn = db + METAGENOME_FILE_SUFFIX;
	string fmidxFn = db + FMINDEX_FILE_SUFFIX;

	ifstream mtgIn;
	ifstream fmidxIn;

	MetaGenome mtg;
	FMIndex fmidx;

	mtgIn.open(mtgFn.c_str(), ios_base::binary);
	if(!mtgIn.is_open()) {
		cerr << "Unable to open '" << mtgFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	fmidxIn.open(fmidxFn.c_str(), ios_base::binary);
	if(!fmidxIn.is_open()) {
		cerr << "Unable to open '" << fmidxFn << "': " << ::strerror(errno) << endl;
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

	infoLog << "Loading FM-index ..." << endl;
	loadProgInfo(fmidxIn);
	if(!fmidxIn.bad())
		fmidx.load(fmidxIn);
	if(fmidxIn.bad()) {
		cerr << "Unable to load reference FM-index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* initiate RNG */
	RNG rng(seed);

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
		return main_SE(mtg, fmidx, fwdI, out, rng, strand, maxMems, bestFrac, maxReport);
	else
		return main_PE(mtg, fmidx, fwdI, revI, out, rng, strand, maxMems, bestFrac, maxReport);
}

int main_SE(const MetaGenome& mtg, const FMIndex& fmidx,
		SeqIO& seqI, SAMfile& out, RNG& rng, int strand,
		uint32_t maxMems, double bestFrac, uint32_t maxReport) {
	infoLog << "Aligning input reads" << endl;
	const DNAseq& target = mtg.getSeq(); // target is always the entire metagenome

	/* search MEMS for each read */
#pragma omp parallel
	{
#pragma omp single
		{
			while(seqI.hasNext()) {
				PrimarySeq read = seqI.nextSeq();
#pragma omp task firstprivate(read)
				{
					MEMS mems = MEMS::sampleMEMS(&read, &fmidx, rng, strand);
					mems.findLocs(); // only find locs for best MEMS
					const MEM::STRAND readStrand = mems.getStrand();
					if(readStrand == MEM::REV)
						read.revcom();

					/* get all candidate MatchPairs */
					const string& id = read.getName();
					const DNAseq& query = read.getSeq();
					const QualStr& qual = read.getQual();
					Alignment::SeedMatchList seedMatches = Alignment::getSeedMatchList(mtg, mems, maxMems);
					if(seedMatches.empty()) {
						infoLog << "Unable to find any valid SeedMatches for '" << id << "', ignore" << endl;
					}
					else {
						vector<Alignment> alnList;
						alnList.reserve(seedMatches.size());
						for(Alignment::SeedMatch& seedMatch : seedMatches) {
							/* get candidate region of this seedMatch */
							int32_t tid = seedMatch.getTId();
							const string& chrName = mtg.getChromName(tid);

							uint64_t tStart = seedMatch.getStart() - seedMatch.getFrom() * (1 + Alignment::MAX_INDEL_RATE);
							uint64_t tEnd = seedMatch.getEnd() + (query.length() - seedMatch.getTo()) * (1 + Alignment::MAX_INDEL_RATE);
							if(tStart < mtg.getChromStart(tid)) // searhStart too far
								tStart = mtg.getChromStart(tid);
							if(tEnd > mtg.getChromEnd(tid)) // searhEnd too far
								tEnd = mtg.getChromEnd(tid);
//							cerr << "id: " << id << " tid: " << tid << " tStart: " << tStart << " tEnd: " << tEnd << endl;
//							cerr << "query:  " << query << endl << "target: " << target.substr(tStart, tEnd - tStart) << endl;
							/* add a new alignment */
							alnList.push_back(Alignment(&query, &target, &qual, &id, tid,
									0, query.length(), tStart, tEnd, &ss, (readStrand == MEM::FWD ? 0 : BAM_FREVERSE)
							).calculateScores(seedMatch).backTrace().clearScores());
						}
						/* find best alignment by alnScore */
						vector<Alignment>::const_iterator bestAln = std::max_element(alnList.cbegin(), alnList.cend(), [] (const Alignment& lhs, const Alignment& rhs) { return lhs.alnScore > rhs.alnScore; });
						/* first-pass filter alignment by alnScore */
						alnList.erase(std::remove_if(alnList.begin(), alnList.end(), [&] (const Alignment& aln) { return aln.alnScore < bestAln->alnScore * bestFrac; }), alnList.end());
						/* evaluate candidate "best-spectrum" alignments */
						for(Alignment& aln : alnList)
							aln.evaluate();
						/* calculate postP/mapQ for candidates */
						Alignment::calcMapQ(alnList);
						/* sort candidates by their log10-liklihood descreasingly (same order as postP with uniform prior */
						std::sort(alnList.begin(), alnList.end(), [] (const Alignment& lhs, const Alignment& rhs) { return lhs.log10P > rhs.log10P; });
						/* export BAM records */
						size_t numReport = maxReport == 0 ? alnList.size() : std::min((size_t) maxReport, alnList.size());
						for(size_t i = 0; i < numReport; ++i) {
							const Alignment& aln = alnList[i];
							/* construct BAM */
							BAM bamAln = aln.exportBAM();
							/* set standard aux tags */
							bamAln.setAux(NUM_REPORTED_ALIGNMENT_TAG, numReport);
							bamAln.setAux(NUM_TOTAL_ALIGNMENT_TAG, alnList.size());
							bamAln.setAux(MISMATCH_POSITION_TAG, aln.getAlnMDTag());
							/* set customized aux tags */
							bamAln.setAux(ALIGNMENT_LOG10LIK_TAG, aln.log10P);
							bamAln.setAux(ALIGNMENT_POSTERIOR_PROB_TAG, aln.postP);
							/* write BAM */
#pragma omp critical(WRITE_BAM)
							out.write(bamAln);
						}
					}
				} /* end task */
			} /* end each read */
		} /* end single, implicit barrier */
	} /* end parallel */
	return out.close();
}

int main_PE(const MetaGenome& mtg, const FMIndex& fmidx,
		SeqIO& fwdI, SeqIO& revI, SAMfile& out, RNG& rng, int strand,
		uint32_t maxMems, double bestFrac, uint32_t maxReport) {
	return EXIT_SUCCESS;
}
