/*
 * BitSeqGGMN.cpp
 *
 *  Created on: Sep 24, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include <cassert>
#include "BitSeqGGMN.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

BitSeqGGMN::BitSeqGGMN(const BitSeqGGMN& other) :
		bstr(other.bstr), factor(other.factor), b(other.b), s(other.s)
{
	size_t nSB = other.numSuperBlocks();
	Rs = new size_t[nSB + 1];
	std::copy(Rs, Rs + nSB, other.Rs);
}

BitSeqGGMN::BitSeqGGMN(const BitStr<data_type>& bstr, size_t factor) : bstr(bstr) {
	if(factor == 0)
		factor = bits(length() - 1);
	this->factor = factor;
	b = sizeof(data_type) * Wb;
	s = b * factor;
	buildRank();
}

void BitSeqGGMN::buildRank() {
	size_t nSB = numSuperBlocks();
	Rs = new size_t[nSB + 1](); /* value initiation */
	for (size_t i = 1; i <= nSB; ++i)
		Rs[i] = Rs[i-1] + buildRank((i-1) * factor, factor);
}

size_t BitSeqGGMN::getBytes() const {
	return bstr.getBytes() + sizeof(size_t) * (numSuperBlocks() + 1) +
			sizeof(factor) + sizeof(b) + sizeof(s) + sizeof(this);
}

size_t BitSeqGGMN::rank1(size_t i) const {
	assert(bstr.getWid() == W);
	size_t r = Rs[i / s]; /* ones before this super-block */
	for(size_t b = i / s * factor; b < i / W; ++b)
		r += popcount(bstr.getValue(b)); /* ones before this data block */
	r += popcount(bstr.get(i / W) & ( (1 << (i & lowest5_mask)) - 1)); /* ones in bits */
	return r;
}

size_t BitSeqGGMN::select1(size_t i) const {
	if(i == 0)
		return -1;

	// binary search over super-blocks
	size_t l = 0;
	size_t r = length() / s;
	size_t mid = (l + r) / 2;
	size_t rankmid = Rs[mid];
	while (l <= r) {
		if (rankmid < i)
			l = mid + 1;
		else
			r = mid - 1;
		mid = (l + r ) / 2;
		rankmid = Rs[mid];
	}
	// sequential search using popcount
	size_t left = mid * factor;
	i -= rankmid;
	size_t j = bstr.getValue(left);
	size_t ones = popcount(j);
	while (ones < i) {
		i -= ones;
		left++;
		if(left > bstr.numValues())
			return length();
		j = bstr.getValue(left);
		ones = popcount(j);
	}

	/* sequential search using popcount over lowest bits */
	left = left * b;
	rankmid = popcount8(j);
	while(rankmid < i) {
		j >>= 8;
		i -= rankmid;
		left += 8;
		rankmid = popcount8(j);
	}

	/* then sequential search bit by bit */
	while(i > 0) {
		if (j & 1)
			i--;
		j >>= 1;
		left++;
	}
	return left - 1;
}

size_t BitSeqGGMN::select0(size_t i) const {
	if(i == 0)
		return -1;
	assert(bstr.getWid() == W);

	/* binary search over first level rank structure */
	size_t l = 0;
	size_t r = length() / s;
	size_t mid = (l + r) / 2;
	size_t rankmid = mid * factor * W - Rs[mid];
	while (l <= r) {
		if (rankmid < i)
			l = mid + 1;
		else
			r = mid - 1;
		mid = (l + r) / 2;
		rankmid = mid * factor * W - Rs[mid];
	}
	/* sequential search using popcount */
	size_t left = mid * factor;
	i -= rankmid;
	size_t j = bstr.getValue(left);
	size_t zeros = W - popcount(j);
	while(zeros < i) {
		i -= zeros;
		left++;
		if(left > bstr.numValues())
			return length();
		j = bstr.getValue(left);
		zeros = W - popcount(j);
	}
	/* sequential search using popcount over lowest bits */
	left = left * b;
	rankmid = 8 - popcount8(j);
	while(rankmid < i) {
		j >>= 8;
		i -= rankmid;
		left += 8;
		rankmid = 8 - popcount8(j);
	}

	/* then sequential search bit by bit */
	while (i > 0) {
		if (j % 2 == 0)
			i--;
		j >>= 1;
		left++;
	}
	left--;
	return left > length() ? length() : left;
}

size_t BitSeqGGMN::selectNext1(size_t start) const {
	assert(bstr.getWid() == W);
	size_t des = start % W;
	uint aux = bstr.getValue(start / W) >> des;

	if (aux) { /* has reminding */
		if(aux& 0xff)
			return start + select_tab[aux & 0xff] - 1;
		else if(aux & 0xff00 > 0)
			return start + 8 + select_tab[(aux >> 8) & 0xff] - 1;
		else if(aux & 0xff0000)
			return start + 16 + select_tab[(aux >> 16) & 0xff] - 1;
		else
			return start + 24 + select_tab[(aux >> 24) & 0xff] - 1;
	}

	/* no reminding */
	for(uint i = start / W + 1; i < bstr.numValues(); ++i) {
		aux = bstr.getValue(i);
		if(aux) {
			uint shift = 0;
			while((aux & (0xff << shift)) == 0)
				shift += Wb;
			return i * W + shift + select_tab[(aux >> shift) & 0xff] - 1;
		}
	}
	return length();
}

size_t BitSeqGGMN::selectPrev1(size_t start) const {
	assert(bstr.getWid() == W);

	size_t i = start >> 5;
	size_t offset = start % W;
	size_t aux = bstr.getValue(i) & (-1u >> (31 - offset));

	if (aux) { /* has remaining */
		uint shift = 0;

		if ((aux&0xFF000000) > 0) return i*W+23+prev_tab[(aux>>24)&0xFF];
		else if ((aux&0xFF0000) > 0) return i*W+15+prev_tab[(aux>>16)&0xFF];
		else if ((aux&0xFF00) > 0) return i*W+7+prev_tab[(aux2>>8)&0xFF];
		else  return i*W+prev_tab[aux&0xFF]-1;
	}
	for (uint k=i-1;;k--) {
		aux=data[k];
		if (aux > 0) {
			if ((aux&0xFF000000) > 0) return k*W+23+prev_tab[(aux>>24)&0xFF];
			else if ((aux&0xFF0000) > 0) return k*W+15+prev_tab[(aux>>16)&0xFF];
			else if ((aux&0xFF00) > 0) return k*W+7+prev_tab[(aux>>8)&0xFF];
			else  return k*W+prev_tab[aux&0xFF]-1;
		}
	}
	return 0;
}

ostream& BitSeqGGMN::save(ostream& out) const {
}

istream& BitSeqGGMN::load(istream& in) {
}

size_t BitSeqGGMN::buildRank(size_t start, size_t len) {
	size_t rank = 0;
	for(size_t i = start; i < bstr.numValues() && i < start + len; ++i)
		rank += popcount32(bstr.getValue(i));
	return rank;
}

} /* namespace libSDS */
} /* namespace EGriceLab */

