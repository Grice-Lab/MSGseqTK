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
 *  A => 1, C => 2, G => 3, T => 4, N => 5 and all invalid symbols as 0
 *  the assembly gaps (Ns) are treated specially
 *  but accept IUPAC ambiguous codes except the gaps
 */
class DNAalphabet {
public:
	/* nested enum and types */
	enum Base { INVALID, A, C, G, T, N };

private:
	/* static fields */
	static const uint8_t* sym2base; /* internal map from symbol to base, static zero initiated by default */
	static const char* base2sym; /* internal map from base to symbol, static zero initiated by default */
	static const char* sym2comp;  /* internal complement map from symbol to symbol, static zero initiated by default */
	static const uint8_t* base2comp;  /* internal complement map from base to base, static zero initiated by default */

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

	/** test whether a base is valid */
	static bool isValid(uint8_t b) {
		return b > 0;
	}

	/** test whether a base is a valid base */
	static bool isBase(uint8_t b) {
		return b != INVALID && b != N;
	}

	/** test whether a base is a gap */
	static bool isGap(uint8_t b) {
		return b == N;
	}

	/** test whether a symbol is valid */
	static bool isValid(char s) {
		return isValid(encode(s));
	}

	/** test whether a base is a valid base */
	static bool isBase(char s) {
		return isBase(encode(s));
	}

	/** test whether a base is a gap */
	static bool isGap(char s) {
		return isGap(encode(s));
	}

	/** get a reverse-completed copy of a dna string */
	static string revcom(const string& str);

private:
	/** initiate the internal sym2base map */
	static uint8_t* initSym2Base();

	/** initiate the internal base2sym map */
	static char* initBase2Sym();

	/** initiate the internal sym2comp map */
	static char* initSym2Comp();

	/** initiate the internal base2comp map */
	static uint8_t* initBase2Comp();

};

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNAALPHABET_H_ */
