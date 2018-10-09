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
		BitStr<uint8_t> bst(Wb, v);
		cout << "bstr(" << (int) v << "): " << bst << endl;
	}

	/* test C-array construction */
	const uint8_t hex_val[] = {
			0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF
	};
	const size_t n_hex = sizeof(hex_val) / sizeof(hex_val[0]);

	BitStr<uint8_t> bst_hex(hex_val, n_hex);
	cerr << "n_hex: " << n_hex << " bst_hex: " << bst_hex << endl;
	for(size_t i = 0; i < n_hex; ++i) {
		fprintf(stdout, "bst_hex.getValue(%d): %ld <-> hex_val[%d]: %ld\n", i, bst_hex.getValue(i), i, hex_val[i]);
		if(bst_hex.getValue(i) != hex_val[i])
			return EXIT_FAILURE;
	}
	/* bitwise test */
	BitStr<uint32_t> bst(W);
	if(bst.getValue(0) != 0)
		return EXIT_FAILURE;
	for(size_t i = 0, v = 1; i < n_hex; ++i, v = v << 1) {
		bst.set(i);
		fprintf(stdout, "bst.getValue(0): %d <-> v: %d\n", bst.getValue(0), v);
		if(bst.getValue(0) != v)
			return EXIT_FAILURE;
		bst.reset(i);
	}

	/* flip test */
	bst.clear();
	for(size_t i = 0, v = 1; i < n_hex; ++i, v = v << 1) {
		bst.flip(i);
		fprintf(stdout, "bst.getValue(0): %d <-> v: %d\n", bst.getValue(0), v);
		if(bst.getValue(0) != v)
			return EXIT_FAILURE;
		bst.flip(i);
	}

	/* logit tests */
	bst.resize(16);
	bst.clear();
	if(bst.any())
		return EXIT_FAILURE;
	cout << "bst.length(): " << bst.length() << endl;
	for(size_t i = 0; i < bst.length(); ++i) {
		bst.set(i);
		fprintf(stdout, "bst.count(): %d i+1: %d\n", bst.count(), i + 1);
		if(bst.count() != i + 1)
			return EXIT_FAILURE;
	}
	fprintf(stdout, "bst.getValue(0): %u\n", bst.getValue(0));
	if(!bst.all())
		return EXIT_FAILURE;

	/* resize and assignment */
	bst.resize(1024);
	cout << "bst.numBits(): " << bst.numBits() << endl;
	if(bst.getValue(0) != UINT16_MAX)
		return EXIT_FAILURE;
	bst.resize(128);
	cout << "bst.numBits(): " << bst.numBits() << endl;
	if(bst.getValue(0) != UINT16_MAX)
		return EXIT_FAILURE;
	bst.resize(4);
	cout << "bst.numBits(): " << bst.numBits() << endl;
	if(bst.getValue(0) != 0xF)
		return EXIT_FAILURE;

	BitStr<uint32_t> bst1 = bst;
	cout << "bst1.numBits(): " << bst1.numBits() << endl;
	if(bst1 != bst)
		return EXIT_FAILURE;
	bst1.clear();
	cout << "bst1: " << bst1 << endl << "bst: " << bst << endl;
	if(!(bst1.none() && bst.any()))
			return EXIT_FAILURE;

	/* copy between different types */
	bst.resize(256);
	BitStr<uint16_t> bstS(bst);
	cout << "bst.getWid(): " << bst.getWid() << " bstS.getWid(): " << bstS.getWid() << endl;
	cout << "bst.numBits(): " << bst.numBits() << " bstS.numBits(): " << bstS.numBits() << endl;
	cout << "bst : " << bst << endl << "bstS: " << bstS << endl;
	if(bst.numBits() != bstS.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst.length(); ++i)
		if(bst.test(i) != bstS.test(i))
			return EXIT_FAILURE;

	BitStr<uint64_t> bstL(bst);
	cout << "bst.getWid(): " << bst.getWid() << " bstL.getWid(): " << bstL.getWid() << endl;
	cout << "bst.numBits(): " << bst.numBits() << " bstL.numBits(): " << bstL.numBits() << endl;
	cout << "bst : " << bst << endl << "bstL: " << bstL << endl;
	if(bst.numBits() != bstL.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst.length(); ++i)
		if(bst.test(i) != bstL.test(i))
			return EXIT_FAILURE;

	/* assign between different types */
	bst1 = BitStr<uint32_t>(256, 0xFF);
	BitStr<uint16_t> bst1S = bst1;
	cout << "bst1.getWid(): " << bst1.getWid() << " bst1S.getWid(): " << bst1S.getWid() << endl;
	cout << "bst1.numBits(): " << bst1.numBits() << " bst1S.numBits(): " << bst1S.numBits() << endl;
	cout << "bst1 : " << bst1 << endl << "bst1S: " << bst1S << endl;
	if(bst1.numBits() != bst1S.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst1.length(); ++i)
		if(bst1.test(i) != bst1S.test(i))
			return EXIT_FAILURE;

	BitStr<uint64_t> bst1L = bst1;
	cout << "bst1.getWid(): " << bst1.getWid() << " bst1L.getWid(): " << bst1L.getWid() << endl;
	cout << "bst1.numBits(): " << bst1.numBits() << " bst1L.numBits(): " << bst1L.numBits() << endl;
	cout << "bst1 : " << bst1 << endl << "bst1L: " << bst1L << endl;
	if(bst1.numBits() != bst1L.numBits())
		return EXIT_FAILURE;
	for(size_t i = 0; i < bst1.length(); ++i)
		if(bst1.test(i) != bst1L.test(i))
			return EXIT_FAILURE;
}
