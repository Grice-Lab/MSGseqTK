/*
 * BitSeq.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#include "BitSeq.h"

namespace EGriceLab {
namespace libSDS {

size_t BitSeq::rank0(const size_t i) const {
	return i + 1 - rank1(i);
}

size_t BitSeq::selectNext1(const size_t i) const
{
	return select1( (i == 0 ? 0 : rank1(i - 1)) + 1);
}

size_t BitSeq::selectNext0(const size_t i) const
{
	return select0( (i == 0 ? 0 : rank0(i - 1)) + 1);
}

size_t BitSeq::selectPrev1(const size_t i) const
{
	size_t v = rank1(i);
	if(v <= 1)
		return -1;
	else
		return select1(v - 1);
}

size_t BitSeq::selectPrev0(const size_t i) const
{
	size_t v = rank0(i);
	if(v <= 1)
		return -1;
	else
		return select0(v - 1);
}

bool BitSeq::access(const size_t i) const
{
	if(i == 0)
		return rank1(i) > 0;
	else
		return rank1(i) > rank1(i - 1);
}

} /* namespace libSDS */
} /* namespace EGriceLab */
