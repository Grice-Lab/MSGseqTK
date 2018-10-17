/*
 * BitSeq.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#include "BitSeq.h"

namespace EGriceLab {
namespace libSDS {

size_t BitSeq::rank0(size_t i) const {
	return i + 1 - rank1(i);
}

size_t BitSeq::selectNext1(size_t start) const
{
	return select1( (start == 0 ? 0 : rank1(start - 1)) + 1);
}

size_t BitSeq::selectNext0(size_t start) const
{
	return select0( (start == 0 ? 0 : rank0(start - 1)) + 1);
}

size_t BitSeq::selectPrev1(size_t start) const
{
	size_t v = rank1(start);
	if(v <= 1)
		return 0;
	else
		return select1(v);
}

size_t BitSeq::selectPrev0(size_t start) const
{
	size_t v = rank0(start);
	if(v <= 1)
		return 0;
	else
		return select0(v);
}

bool BitSeq::access(size_t i) const
{
	if(i == 0)
		return rank1(i) > 0;
	else
		return rank1(i) > rank1(i - 1);
}

ostream& BitSeq::save(ostream& out) const {
	out.write((const char*) &n, sizeof(size_t));
	out.write((const char*) &ones, sizeof(size_t));
	return out;
}

istream& BitSeq::load(istream& in) {
	in.read((char*) &n, sizeof(size_t));
	in.read((char*) &ones, sizeof(size_t));
	return in;
}

} /* namespace libSDS */
} /* namespace EGriceLab */
