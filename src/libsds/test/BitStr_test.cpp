/*
 * BitStr_test.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#include <string>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <cstdio>
#include <climits>
#include <vector>
#include "BitStr.h"

namespace EGriceLab {
namespace libSDS {

} /* namespace libSDS */
} /* namespace EGriceLab */

using namespace std;
using namespace EGriceLab::libSDS;

int main() {
	/* try hex values */
	for(uint8_t v = 0; v < 16; ++v) {
		BitStr8 bst(Wb, v);
		cerr << "bstr(" << (int) v << "): " << bst << endl;
	}

	/* test C-array construction */
	const uint8_t hex_val[] = {
			0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF
	};
	const size_t n_hex = sizeof(hex_val) / sizeof(hex_val[0]);

	BitStr8 bst_hex(hex_val, n_hex);
	cerr << "bst_hex.numValues(): " << bst_hex.numValues() << " n_hex: " << n_hex << " bst_hex: " << bst_hex << endl;
	for(size_t i = 0; i < n_hex; ++i) {
		fprintf(stderr, "bst_hex.getValue(%d): %ld <-> hex_val[%d]: %ld\n", i, bst_hex.getValue(i), i, hex_val[i]);
		if(bst_hex.getValue(i) != hex_val[i])
			return EXIT_FAILURE;
	}
	/* test C-array construction of different type */
	BitStr32 bst_hex32(hex_val, n_hex);
	cerr << "bst_hex32.numValues(): " << bst_hex32.numValues() << " n_hex: " << n_hex << " bst_hex: " << bst_hex << endl;
	for(size_t i = 0; i < n_hex; ++i) {
		fprintf(stderr, "bst_hex.getValue(%d, %d): %ld <-> hex_val[%d]: %ld\n", i, 8, bst_hex.getValue(i, 8), i, hex_val[i]);
		if(bst_hex.getValue(i, 8) != hex_val[i])
			return EXIT_FAILURE;
	}

	/* test std::string construction */
	const string hex_str = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF };
	bst_hex = BitStr8(hex_str);
	cerr << "bst_hex.numValues(): " << bst_hex.numValues() << " n_hex: " << n_hex << " bst_hex: " << bst_hex << endl;
	for(size_t i = 0; i < n_hex; ++i) {
		fprintf(stderr, "bst_hex.getValue(%d): %ld <-> hex_str[%d]: %ld\n", i, bst_hex.getValue(i), i, hex_str[i]);
		if(bst_hex.getValue(i) != hex_str[i])
			return EXIT_FAILURE;
	}

	/* bitwise test */
	BitStr32 bst(W);
	if(bst.getValue(0) != 0)
		return EXIT_FAILURE;
	for(size_t i = 0, v = 1; i < n_hex; ++i, v = v << 1) {
		bst.set(i);
		fprintf(stderr, "bst.getValue(0): %d <-> v: %d\n", bst.getValue(0), v);
		if(bst.getValue(0) != v)
			return EXIT_FAILURE;
		bst.reset(i);
	}

	/* flip test */
	bst.clear();
	for(size_t i = 0, v = 1; i < n_hex; ++i, v = v << 1) {
		bst.flip(i);
		fprintf(stderr, "bst.getValue(0): %d <-> v: %d\n", bst.getValue(0), v);
		if(bst.getValue(0) != v)
			return EXIT_FAILURE;
		bst.flip(i);
	}

	/* logit tests */
	bst.resize(16);
	bst.clear();
	if(bst.any())
		return EXIT_FAILURE;
	cerr << "bst.length(): " << bst.length() << endl;
	for(size_t i = 0; i < bst.length(); ++i) {
		bst.set(i);
		fprintf(stderr, "bst.count(): %d i+1: %d\n", bst.count(), i + 1);
		if(bst.count() != i + 1)
			return EXIT_FAILURE;
	}
	fprintf(stderr, "bst.getValue(0): %u\n", bst.getValue(0));
	if(!bst.all())
		return EXIT_FAILURE;

	/* value setting tests */
	BitStr32 bst1(512);
	for(size_t i = 0, len = 1; i * bst1.getWid() < bst1.length() && len < 16; ++len) {
		uint32_t val = len;
		fprintf(stderr, "bst1.setValue(%d, %d, %d)\n", i, len, val);
		bst1.setValue(i, len, val);
		uint32_t valN = bst1.getValue(i, len);
		fprintf(stderr, "bst1.getValue(%d, %d): %d\n", i, len, valN);
		if(valN != val)
			return EXIT_FAILURE;
		i += len;
	}

	/* bit setting tests */
	bst1.clear();
	for(size_t i = 0, len = 1; i < bst1.length() && len < 32; ++len) {
		uint32_t val = len;
		fprintf(stderr, "bst1.set(%d, %d, %d)\n", i, len, val);
		bst1.set(i, len, val);
		uint32_t valN = bst1.get(i, len);
		fprintf(stderr, "bst1.get(%d, %d): %d\n", i, len, valN);
		if(valN != val)
			return EXIT_FAILURE;
		i += len;
	}

	/* resize and assignment */
	bst.resize(1024);
	cerr << "bst.numBits(): " << bst.numBits() << endl;
	if(bst.getValue(0) != UINT16_MAX)
		return EXIT_FAILURE;
	bst.resize(128);
	cerr << "bst.numBits(): " << bst.numBits() << endl;
	if(bst.getValue(0) != UINT16_MAX)
		return EXIT_FAILURE;
	bst.resize(4);
	cerr << "bst.numBits(): " << bst.numBits() << endl;
	if(bst.getValue(0) != 0xF)
		return EXIT_FAILURE;

	bst1 = bst;
	cerr << "bst1.numBits(): " << bst1.numBits() << endl;
	if(bst1 != bst)
		return EXIT_FAILURE;
	bst1.clear();
	cerr << "bst1: " << bst1 << endl << "bst: " << bst << endl;
	if(!(bst1.none() && bst.any()))
			return EXIT_FAILURE;

	/* copy between different types */
	bst.resize(256);
	BitStr<uint16_t> bstS(bst);
	cerr << "bst.getWid(): " << bst.getWid() << " bstS.getWid(): " << bstS.getWid() << endl;
	cerr << "bst.numBits(): " << bst.numBits() << " bstS.numBits(): " << bstS.numBits() << endl;
	cerr << "bst : " << bst << endl << "bstS: " << bstS << endl;
	if(bst.numBits() != bstS.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst.length(); ++i)
		if(bst.test(i) != bstS.test(i))
			return EXIT_FAILURE;

	BitStr64 bstL(bst);
	cerr << "bst.getWid(): " << bst.getWid() << " bstL.getWid(): " << bstL.getWid() << endl;
	cerr << "bst.numBits(): " << bst.numBits() << " bstL.numBits(): " << bstL.numBits() << endl;
	cerr << "bst : " << bst << endl << "bstL: " << bstL << endl;
	if(bst.numBits() != bstL.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst.length(); ++i)
		if(bst.test(i) != bstL.test(i))
			return EXIT_FAILURE;

	/* assign between different types */
	bst1 = BitStr<uint32_t>(256, 0xFF);
	BitStr16 bst1S = bst1;
	cerr << "bst1.getWid(): " << bst1.getWid() << " bst1S.getWid(): " << bst1S.getWid() << endl;
	cerr << "bst1.numBits(): " << bst1.numBits() << " bst1S.numBits(): " << bst1S.numBits() << endl;
	cerr << "bst1 : " << bst1 << endl << "bst1S: " << bst1S << endl;
	if(bst1.numBits() != bst1S.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst1.length(); ++i)
		if(bst1.test(i) != bst1S.test(i))
			return EXIT_FAILURE;

	BitStr64 bst1L = bst1;
	cerr << "bst1.getWid(): " << bst1.getWid() << " bst1L.getWid(): " << bst1L.getWid() << endl;
	cerr << "bst1.numBits(): " << bst1.numBits() << " bst1L.numBits(): " << bst1L.numBits() << endl;
	cerr << "bst1 : " << bst1 << endl << "bst1L: " << bst1L << endl;
	if(bst1.numBits() != bst1L.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst1.length(); ++i)
		if(bst1.test(i) != bst1L.test(i))
			return EXIT_FAILURE;
}
