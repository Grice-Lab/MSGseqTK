/*
 * libsdsBitBasic.h
 *  basic non-class functions for bit operations
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#ifndef LIBSDSBITBASIC_H_
#define LIBSDSBITBASIC_H_

#include <cstddef>
#include "libsdsConst.h"

namespace EGriceLab {
namespace libSDS {

/** bits needed to represent a number between 0 and n */
inline uint bits(size_t n) {
	uint b = 0;
	for(; n != 0; n >>= 1)
		b++;
	return b;
}

/** Counts the number of 1s in x */
inline uint popcount(const int x) {
	return __popcount_tab[(x >>  0) & 0xff]  + __popcount_tab[(x >>  8) & 0xff]
															  + __popcount_tab[(x >> 16) & 0xff] + __popcount_tab[(x >> 24) & 0xff];
}

/** Counts the number of 1s in the first 16 bits of x */
inline uint popcount16(const int x) {
	return __popcount_tab[x & 0xff]  + __popcount_tab[(x >>  8) & 0xff];
}

/** Counts the number of 1s in the first 8 bits of x */
inline uint popcount8(const int x) {
	return __popcount_tab[x & 0xff];
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* LIBSDSBITBASIC_H_ */
