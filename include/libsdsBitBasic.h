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
inline unsigned int bits(size_t n) {
	size_t b = 0;
	for(; n != 0; n >>= 1)
		b++;
	return b;
}

/** counts the number of 1s in x */
template<typename uIntType>
inline unsigned int popcount(uIntType x) {
	unsigned int ones = 0;
	for(size_t i = 0; i < sizeof(uIntType); ++i)
		ones += __popcount_tab[(x >> Wb * i) & 0xff];
	return ones;
}

/** Counts the number of 1s in the first 64 bits of x */
inline unsigned int popcount64(uint64_t x) {
	return __popcount_tab[(x >> 0) & 0xff] + __popcount_tab[(x >> 8) & 0xff]
         + __popcount_tab[(x >> 16) & 0xff] + __popcount_tab[(x >> 24) & 0xff]
		 + __popcount_tab[(x >> 32) & 0xff] + __popcount_tab[(x >> 40) & 0xff]
		 + __popcount_tab[(x >> 48) & 0xff] + __popcount_tab[(x >> 56) & 0xff];
}

/** Counts the number of 1s in the first 32 bits of x */
inline unsigned int popcount32(uint32_t x) {
	return __popcount_tab[(x >> 0) & 0xff] + __popcount_tab[(x >> 8) & 0xff]
         + __popcount_tab[(x >> 16) & 0xff] + __popcount_tab[(x >> 24) & 0xff];
}

/** Counts the number of 1s in the first 16 bits of x */
inline unsigned int popcount16(uint16_t x) {
	return __popcount_tab[x & 0xff]  + __popcount_tab[(x >> 8) & 0xff];
}

/** Counts the number of 1s in the first 8 bits of x */
inline unsigned int popcount8(uint8_t x) {
	return __popcount_tab[x & 0xff];
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* LIBSDSBITBASIC_H_ */
