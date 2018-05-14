/*
 * DNAalphabet.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "DNAalphabet.h"

namespace EGriceLab {
namespace MSGseqClean {

const int8_t* DNAalphabet::sym2base = initSym2Base();
const char* DNAalphabet::base2sym = initBase2Sym();
const char* DNAalphabet::sym2comp = initSym2Comp();
const int8_t* DNAalphabet::base2comp = initBase2Comp();

int8_t* DNAalphabet::initSym2Base() {
	/* sym2base will be zero initiated by default */
	static int8_t sym2base[CHAR_MAX + 1];
	std::fill_n(sym2base, CHAR_MAX + 1, -1);
	/* basic symbols */
	sym2base['A'] = sym2base['a'] = A;
	sym2base['C'] = sym2base['c'] = C;
	sym2base['G'] = sym2base['g'] = G;
	sym2base['T'] = sym2base['t'] = T;
	sym2base['N'] = sym2base['n'] = N;
	/* IUPAC extensions are resolved arbitrary */
	sym2base['U'] = sym2base['u'] = sym2base['T'];
	sym2base['R'] = sym2base['r'] = sym2base['A']; /* R -> A/G */
	sym2base['Y'] = sym2base['y'] = sym2base['C']; /* Y -> C/T */
	sym2base['S'] = sym2base['s'] = sym2base['G']; /* S -> C/G */
	sym2base['W'] = sym2base['w'] = sym2base['T']; /* W -> A/T */
	sym2base['K'] = sym2base['k'] = sym2base['G']; /* K -> G/T */
	sym2base['M'] = sym2base['m'] = sym2base['C']; /* M -> A/C */
	sym2base['B'] = sym2base['b'] = sym2base['T']; /* B -> C/G/T */
	sym2base['D'] = sym2base['d'] = sym2base['A']; /* D -> A/G/T */
	sym2base['H'] = sym2base['h'] = sym2base['C']; /* H -> A/C/T */
	sym2base['V'] = sym2base['v'] = sym2base['G']; /* V -> A/C/G */
	return sym2base;
}

char* DNAalphabet::initBase2Sym() {
	/* base2sym will be zero initiated by default */
	static char base2sym[INT8_MAX + 1];
	base2sym[A] = 'A';
	base2sym[C] = 'C';
	base2sym[G] = 'G';
	base2sym[T] = 'T';
	base2sym[N] = 'N';
	return base2sym;
}

char* DNAalphabet::initSym2Comp() {
	static char sym2comp[CHAR_MAX]; /* sym2comp  will be zero initiated */
	/* basic symbols, upper case */
	sym2comp['N'] = 'N';
	sym2comp['A'] = 'T';
	sym2comp['T'] = 'A';
	sym2comp['C'] = 'G';
	sym2comp['G'] = 'C';
	/* basic symbols, lower case */
	sym2comp['n'] = 'n';
	sym2comp['a'] = 't';
	sym2comp['t'] = 'a';
	sym2comp['c'] = 'g';
	sym2comp['g'] = 'c';
	/* IUPAC extensions, upper case */
	sym2comp['U'] = 'A';
	sym2comp['Y'] = 'R';
	sym2comp['R'] = 'Y';
	sym2comp['S'] = 'S';
	sym2comp['W'] = 'W';
	sym2comp['K'] = 'M';
	sym2comp['M'] = 'K';
	sym2comp['B'] = 'V';
	sym2comp['V'] = 'B';
	sym2comp['D'] = 'H';
	sym2comp['H'] = 'D';
	/* IUPAC extensions, lower case */
	sym2comp['u'] = 'a';
	sym2comp['y'] = 'r';
	sym2comp['r'] = 'y';
	sym2comp['s'] = 's';
	sym2comp['w'] = 'w';
	sym2comp['k'] = 'm';
	sym2comp['m'] = 'k';
	sym2comp['b'] = 'v';
	sym2comp['v'] = 'b';
	sym2comp['d'] = 'h';
	sym2comp['h'] = 'd';
	return sym2comp;
}

int8_t* DNAalphabet::initBase2Comp() {
	static int8_t base2comp[INT8_MAX];
	base2comp[A] = T;
	base2comp[C] = G;
	base2comp[G] = C;
	base2comp[T] = A;
	base2comp[N] = N;
	return base2comp;
}

string DNAalphabet::revcom(const string& str) {
	string rcStr;
	rcStr.reserve(str.length());
	for(string::const_reverse_iterator s = str.rbegin(); s != str.rend(); ++s)
		rcStr.push_back(complement(*s));
	return rcStr;
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */


