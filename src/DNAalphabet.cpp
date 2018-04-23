/*
 * DNAalphabet.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "DNAalphabet.h"

namespace EGriceLab {
namespace MSGSeqClean {

void DNAalphabet::initSym2Base(uint8_t* sym2base) {
	/* basic symbols */
	sym2base['N'] = sym2base['n'] = 0;
	sym2base['A'] = sym2base['a'] = 1;
	sym2base['C'] = sym2base['c'] = 2;
	sym2base['G'] = sym2base['g'] = 3;
	sym2base['T'] = sym2base['t'] = 4;
	/* IUPAC extensions */
	sym2base['U'] = sym2base['u'] = sym2base['T'];
	sym2base['R'] = sym2base['r'] = sym2base['A']; /* R -> A/G */
	sym2base['Y'] = sym2base['y'] = sym2base['C']; /* Y -> C/T */
	sym2base['S'] = sym2base['s'] = sym2base['C']; /* S -> C/G */
	sym2base['W'] = sym2base['w'] = sym2base['A']; /* W -> A/T */
	sym2base['K'] = sym2base['k'] = sym2base['G']; /* K -> G/T */
	sym2base['M'] = sym2base['m'] = sym2base['A']; /* M -> A/C */
	sym2base['B'] = sym2base['b'] = sym2base['C']; /* B -> C/G/T */
	sym2base['D'] = sym2base['d'] = sym2base['A']; /* D -> A/G/T */
	sym2base['H'] = sym2base['h'] = sym2base['A']; /* H -> A/C/T */
	sym2base['V'] = sym2base['v'] = sym2base['A']; /* V -> A/C/G */
}

void DNAalphabet::initBase2Sym(char* base2sym) {
	std::fill_n(base2sym, UINT8_MAX, 'N'); /* everything except 1..4 will be an N' */
	base2sym[1] = 'A';
	base2sym[2] = 'C';
	base2sym[3] = 'G';
	base2sym[4] = 'T';
}

void DNAalphabet::initSym2Comp(char* sym2comp) {
	/* basic symbols, upper case */
	sym2comp['N'] = 'N';
	sym2comp['A'] = 'T';
	sym2comp['T'] = 'A';
	sym2comp['C'] = 'G';
	sym2comp['G'] = 'C';
	/* basic symbols, lower case */
	sym2comp['n'] = 'n';
	sym2comp['a'] = 'a';
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
}

void DNAalphabet::initBase2Comp(uint8_t* base2comp) {
	base2comp[1] = 4; /* A -> T */
	base2comp[2] = 3; /* C -> G */
	base2comp[3] = 2; /* G -> C */
	base2comp[4] = 1; /* T -> A */
}

} /* namespace MSGSeqClean */
} /* namespace EGriceLab */


