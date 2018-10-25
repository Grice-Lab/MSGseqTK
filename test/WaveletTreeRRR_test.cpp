/*
 * SeqRRR_test.cpp
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cerrno>
#include "WaveletTreeRRR.h"

namespace EGriceLab {
namespace libSDS {

} /* namespace libSDS */
} /* namespace EGriceLab */

using namespace std;
using namespace EGriceLab::libSDS;

int main() {
	const size_t N = 20;
	const string seq = "ACGTNCGTANGTACNTACGN";
	const basic_string<uint8_t> dna = { 0, 1, 2, 3, 4, 1, 2, 3, 0, 4, 2, 3, 0, 1, 4, 3, 0, 1, 2, 4 };
	assert(seq.length() == N);
	assert(dna.length() == N);

	WaveletTreeRRR seqRRR(seq);
	cerr << "seqRRR built" << endl;

	WaveletTreeRRR dnaRRR(dna);
	cerr << "dnaRRR built" << endl;

	/* access tests */
	for(size_t i = 0; i < N; ++i) {
		char c = seqRRR.access(i);
		fprintf(stderr, "seqRRR.access(%d): %c seq[%d]: %c\n", i, c, i, seq[i]);
		if(c != seq[i])
			return EXIT_FAILURE;
	}
	for(size_t i = 0; i < N; ++i) {
		uint8_t c = dnaRRR.access(i);
		fprintf(stderr, "dnaRRR.access(%d): %d seq[%d]: %d\n", i, c, i, dna[i]);
		if(dnaRRR.access(i) != dna[i])
			return EXIT_FAILURE;
	}
	/* rank tests */
	for(char c : string("ACGTN")) {
		for(size_t i = 0; i < N; ++i) {
			size_t r = seqRRR.rank(c, i);
			fprintf(stderr, "seqRRR.rank(%c, %d): %d\n", c, i, r);
			if(r != std::count(seq.begin(), seq.begin() + i + 1, c))
				return EXIT_FAILURE;
		}
	}
	for(uint8_t c : basic_string<uint8_t>{0, 1, 2, 3, 4}) {
		for(size_t i = 0; i < N; ++i) {
			size_t r = dnaRRR.rank(c, i);
			fprintf(stderr, "dnaRRR.rank(%d, %d): %d\n", c, i, r);
			if(r != std::count(dna.begin(), dna.begin() + i + 1, c))
				return EXIT_FAILURE;
		}
	}
	/* select tests */
	for(char c : string("ACGTN")) {
		for(size_t r = 1; r <= N; ++r) {
			size_t i = seqRRR.select(c, r);
			fprintf(stderr, "seqRRR.select(%c, %d): %d\n", c, r, i);
			if(!(i == N || r == seqRRR.rank(c, i)))
				return EXIT_FAILURE;
		}
	}
	for(char c : basic_string<uint8_t>{0, 1, 2, 3, 4}) {
		for(size_t r = 1; r <= N; ++r) {
			size_t i = dnaRRR.select(c, r);
			fprintf(stderr, "dnaRRR.select(%d, %d): %d\n", c, r, i);
			if(!(i == N || r == dnaRRR.rank(c, i)))
				return EXIT_FAILURE;
		}
	}
	/* copy test */
	WaveletTreeRRR seqRRR2 = seqRRR;
	if(seqRRR2 != seqRRR) {
		cerr << "failed to copy seqRRR" << endl;
		return EXIT_FAILURE;
	}
	WaveletTreeRRR dnaRRR2 = dnaRRR;
	if(dnaRRR2 != dnaRRR) {
		cerr << "failed to copy dnaRRR" << endl;
		return EXIT_FAILURE;
	}
	/* IO tests */
	ostringstream out;
	seqRRR.save(out);
	if(out.bad()) {
		cerr << "failed to save seqRRR: '" << ::strerror(errno) << " '" << endl;
		return EXIT_FAILURE;
	}
	istringstream in(out.str());
	seqRRR2.load(in);
	if(in.bad()) {
		cerr << "failed to load seqRRR: '" << ::strerror(errno) << " '" << endl;
		return EXIT_FAILURE;
	}
	if(seqRRR2 != seqRRR) {
		cerr << "loaded seqRRR doesn't match saved object" << endl;
		return EXIT_FAILURE;
	}
	ostringstream out2;
	dnaRRR.save(out2);
	if(out2.bad()) {
		cerr << "failed to save dnaRRR: '" << ::strerror(errno) << " '" << endl;
		return EXIT_FAILURE;
	}
	istringstream in2(out2.str());
	dnaRRR2.load(in2);
	if(in2.bad()) {
		cerr << "failed to load dnaRRR: '" << ::strerror(errno) << " '" << endl;
		return EXIT_FAILURE;
	}
	if(dnaRRR2 != dnaRRR) {
		cerr << "loaded dnaRRR doesn't match saved object" << endl;
		return EXIT_FAILURE;
	}
}
