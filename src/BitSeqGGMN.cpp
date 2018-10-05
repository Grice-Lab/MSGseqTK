/*
 * BitSeqGGMN.cpp
 *
 *  Created on: Sep 24, 2018
 *      Author: zhengqi
 */

#include "BitSeqGGMN.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

BitSeqGGMN::BitSeqGGMN(const BitStr<uint>& bstr, size_t factor) : bstr(bstr) {
	if(factor == 0)
		factor = bits(length() - 1);
	this->factor = factor;
	b = sizeof(uint);
	s = b * factor;
	BuildRank();
}

BitSeqGGMN::BitSeqGGMN(const BitStr<class oIntType>& bst, size_t factor) {
	if(factor == 0)
		factor = bits(length() - 1);
	this->factor = factor;
	/* copy BitStr data */
}

} /* namespace libSDS */
} /* namespace EGriceLab */

