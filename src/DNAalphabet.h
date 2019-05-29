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
#include <cstdint>
#include <array>

namespace EGriceLab {
namespace MSGseqTK {
using std::string;
using std::basic_string;
using std::array;

typedef uint8_t nt16_t;

/**
 *  This DNA alphabet class use the same nt16 (4bits) encoding as in the htslib package except for Ns,
 *  which encodes A, C, G, T as 1/2/4/8 bits, and IUPAC codes as the mixture of these bits,
 *  and encode N as genomic gap 0 (null char),
 *  and uses the null terminal ('\0') as genomic gaps (between chromsomes)
 *  It provide static methods for encoding and decoding characters to/from small integers
 */
class DNAalphabet {
public:
	typedef array<nt16_t, 256> base_map; /* base map whose value type is base */
	typedef array<char, 256> sym_map; /* symbol map whose value type is symbol */
	typedef array<int, 256> int_map; /* base/symbol map whose value type is (small) int */
	/* nested enum and types */
	enum Base
	{   A = 1, C = 2, G = 4, T = 8,
		U = T,
		R = A | G, // R is A or G
		Y = C | T, // Y is C or T
		S = G | C, // S is G or C
		W = A | T, // W is A or T
		K = G | T, // K is G or T
		M = A | C, // M is A or C
		B = C | G | T, // B is C or G or T
		D = A | G | T, // D is A or G or T
		H = A | C | T, // H is A or C or T
		V = A | C | G, // V is A or C or G
		N = A | C | G | T // N is any base
	};
	static const nt16_t GAP_BASE = 0;
	static const nt16_t NT16_MIN = 1;
	static const nt16_t NT16_MAX = 0xF;
	static const nt16_t SIZE = NT16_MAX + 1;
//	static const nt16_t BASIC_SIZE = 4; /* A, C, G, T */
	static const char GAP_SYM = '$';
	static const uint32_t NT16_LOWER_MASK = 0xf;
	static const uint32_t NT16_UPPER_MASK = 0xf << 4;

private:
	/* static fields */
	static const base_map sym2base; /* internal symbo->base map, almost the same as the nt16 encoding except Ns */
	static const sym_map base2sym;  /* internal base->sym map, almost the same as nt16 decoding except Ns */
	static const sym_map sym2comp;  /* internal complement map from symbol to symbol, static zero initiated by default */
	static const base_map base2comp;  /* internal complement map from base to base, static zero initiated by default */
	static const int_map base2int; /* internal base->int(0,1,2,3) map, with default value 4 */
	static const int_map sym2int;   /* internal sym->int(0,1,2,3) map, with default value 4 */
	static const base_map base2basic; /* internal base->base map to map IUPAC extended bases to nearest basic bases */

	/* static methods */
public:
	/** encode a DNA symbol to a base value */
	static nt16_t encode(char s) {
		return sym2base[s];
	}

	/** decode a DNA base to from a char symbol */
	static char decode(nt16_t b) {
		return base2sym[b];
	}

	/** get complement symbol of a given symbol */
	static char complement(char s) {
		return sym2comp[s];
	}

	/** get complement base of a given base */
	static nt16_t complement(nt16_t b) {
		return base2comp[b];
	}

	/** map a base to basic base */
	static nt16_t toBasic(nt16_t b) {
		return base2basic[b];
	}

	/** map a symbol to basic symbol */
	static char toBasic(char s) {
		return DNAalphabet::decode(base2basic[DNAalphabet::encode(s)]);
	}

	/** test whether a base is valid ( a gap or base ) */
	static bool isValid(nt16_t b) {
		return b <= NT16_MAX;
	}

	/** test whether a base is a gap */
	static bool isGap(nt16_t b) {
		return b == GAP_BASE;
	}

	/** test whether a base is a valid base */
	static bool isBase(nt16_t b) {
		return NT16_MIN <= b && b <= NT16_MAX;
	}

	/** test whether a symbol is valid */
	static bool isValid(char s) {
		return isValid(encode(s));
	}

	/** test whether a symbol is a base */
	static bool isBase(char s) {
		return isBase(encode(s));
	}

	/** test whether a symbol is a gap base */
	static bool isGap(char s) {
		return s == GAP_SYM;
	}

	/** test whether a base is ambiguous */
	static bool isAmbiguous(nt16_t b) {
		return !isBasic(b);
	}

	/** test whether a symbol is ambiguous */
	static bool isAmbiguous(char s) {
		return !isBasic(s);
	}

	/** test whether a base is basic base */
	static bool isBasic(nt16_t b) {
		switch(b) {
		case A: case C: case G: case T: case N:
			return true;
		default:
			return false;
		}
	}

	/** test whether a symbol is basic base */
	static bool isBasic(char s) {
		return isBasic(encode(s));
	}

	/** decode a base to small int */
	static int toInt(nt16_t b) {
		return base2int[b];
	}

	/** decode a symbol to small int */
	static int toInt(char s) {
		return sym2int[s];
	}

	/** get a reverse-completed copy of a dna string */
	static string revcom(const string& str);
};

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */

#endif /* SRC_DNAALPHABET_H_ */
