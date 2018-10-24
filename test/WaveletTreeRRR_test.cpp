/*
 * SeqRRR_test.cpp
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#include <cassert>
#include "../include/WaveletTreeRRR.h"

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

//	/* access tests */
//	for(size_t i = 0; i < N; ++i) {
//		char ch = seqRRR.access(i);
//		fprintf(stderr, "seqRRR.access(%d): %c seq[%d]: %c\n", i, ch, i, seq[i]);
////		if(ch != seq[i])
////			return EXIT_FAILURE;
//	}

	for(size_t i = 0; i < N; ++i) {
		uint8_t b = dnaRRR.access(i);
		fprintf(stderr, "dnaRRR.access(%d): %d seq[%d]: %d\n", i, b, i, dna[i]);
//		if(dnaRRR.access(i) != dna[i])
//			return EXIT_FAILURE;
	}
}
