/*
 * BitSeqRRR_test.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cerrno>
#include "BitSeqRRR.h"

namespace EGriceLab {
namespace libSDS {

} /* namespace libSDS */
} /* namespace EGriceLab */

using namespace std;
using namespace EGriceLab::libSDS;

int main() {
	/* construct a small BitStr */
	size_t N = 256;
	BitStr<uint32_t> bstr(N);
	/* set every 4 bits on */
	for(size_t i = 0; i < N; i += 4)
		bstr.set(i);
	cerr << "bstr: " << bstr << endl;
	/* construct a BitSeqRRR */
	BitSeqRRR bseq(bstr);

	/* basic tests */
	cerr << "bseq.length(): " << bseq.length() << " bseq.numOnes(): " << bseq.numOnes() << endl;
	if(!(bseq.length() == N && bseq.numOnes() == N / 4))
		return EXIT_FAILURE;
	/* access tests */
	for(size_t i = 0; i < N; ++i) {
		fprintf(stderr, "bseq.access(%d): %d\n", i, bseq.access(i));
		if(bseq.access(i) != i % 4 == 0)
			return EXIT_FAILURE;
	}
	/* rank tests */
	for(size_t i = 0; i < N; ++i) {
		fprintf(stderr, "bseq.rank1(%d): %d\n", i, bseq.rank1(i));
		if(bseq.rank1(i) != i / 4 + 1)
			return EXIT_FAILURE;
		fprintf(stderr, "bseq.rank0(%d): %d\n", i, bseq.rank0(i));
		if(bseq.rank0(i) != i - i / 4)
			return EXIT_FAILURE;
	}
	/* select tests */
	for(size_t r = 1; r <= bseq.numOnes(); ++r) {
		fprintf(stderr, "bseq.select1(%d): %d\n", r, bseq.select1(r));
		if(bseq.select1(r) != (r - 1) * 4)
			return EXIT_FAILURE;
	}
	for(size_t r = 1; r <= bseq.numZeros(); ++r) {
		fprintf(stderr, "bseq.select0(%d): %d\n", r, bseq.select0(r));
		if(bseq.select0(r) != r + (r - 1) / 3)
			return EXIT_FAILURE;
	}
	/* forward search tests */
	for(size_t i = 0; i < N; ++i) {
		fprintf(stderr, "bseq.selectNext1(%d): %d\n", i, bseq.selectNext1(i));
		if(bseq.selectNext1(i) != (i + 4 - 1) / 4 * 4)
			return EXIT_FAILURE;
	}
	/* reverse search tests */
	for(size_t i = 0; i < N; ++i) {
		fprintf(stderr, "bseq.selectPrev1(%d): %d\n", i, bseq.selectPrev1(i));
		if(bseq.selectPrev1(i) != i / 4 * 4)
			return EXIT_FAILURE;
	}
	/* access & rank test */
	for(size_t i = 0, r = 0; i < N; ++i) {
		fprintf(stderr, "bseq.access(%d, %d): ", i, r);
		bool bit = bseq.access(i, r);
		fprintf(stderr, "%d\n", bit);
		fprintf(stderr, "r: %d\n", r);
		if(bit != i % 4 == 0)
			return EXIT_FAILURE;
		if(r != (bit ? bseq.rank1(i) : i - bseq.rank1(i) + 1))
			return EXIT_FAILURE;
	}

	/* copy test */
	BitSeqRRR bseqN = bseq;
	if(bseqN != bseq) {
		cerr << "failed to copy bseq" << endl;
		return EXIT_FAILURE;
	}

	/* IO tests */
	ostringstream out;
	bseq.save(out);
	if(out.bad()) {
		cerr << "failed to save bseq: '" << ::strerror(errno) << " '" << endl;
		return EXIT_FAILURE;
	}

	istringstream in(out.str());
	BitSeqRRR bseq1;
	bseq1.load(in);
	if(in.bad()) {
		cerr << "failed to load bseq: '" << ::strerror(errno) << " '" << endl;
		return EXIT_FAILURE;
	}
	if(bseq1 != bseq) {
		cerr << "loaded BitSeqRRR doesn't match saved object" << endl;
		return EXIT_FAILURE;
	}
}
