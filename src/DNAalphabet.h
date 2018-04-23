/*
 * DNAalphabet.h
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */

#ifndef SRC_DNAALPHABET_H_
#define SRC_DNAALPHABET_H_

#include <string>
#include <climits>
#include <cstdint> // C++11

namespace EGriceLab {
namespace MSGSeqClean {
using std::string;

/**
 *  A DNA alphabet class provide static methods for encoding and decoding characters to/from small integers
 *  A => 1, C => 2, G => 3, T => 4 and N and invalid symbols as 0,
 *  but accept IUPAC ambiguous codes except the gaps
 */
class DNAalphabet {
	/* static methods */
public:
	/** encode a DNA symbol to a base value */
	static uint8_t encode(char s) {
		return sym2base[s];
	}

	/** decode a DNA base to a symbol */
	static char decode(uint8_t b) {
		return base2sym[b];
	}

	/** get complement symbol of a given symbol */
	static char complement(char s) {
		return sym2comp[s];
	}

	/** get complement base of a given base */
	static uint8_t complement(uint8_t b) {
		return base2comp[b];
	}

private:
	/** initiate the internal sym2base map */
	static void initSym2Base(uint8_t* sym2base);

	/** initiate the internal base2sym map */
	static void initBase2Sym(char* base2sym);

	/** initiate the internal sym2comp map */
	static void initSym2Comp(char* sym2comp);

	/** initiate the internal base2comp map */
	static void initBase2Comp(uint8_t* base2comp);

	/* static fields */
private:
	static uint8_t sym2base[CHAR_MAX]; /* internal map from symbol to base, static zero initiated by default */
	static char base2sym[UINT8_MAX]; /* internal map from base to symbol, static zero initiated by default */
	static char sym2comp[CHAR_MAX];  /* internal complement map from symbol to symbol, static zero initiated by default */
	static uint8_t base2comp[UINT8_MAX];  /* internal complement map from base to base, static zero initiated by default */
};

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNAALPHABET_H_ */
