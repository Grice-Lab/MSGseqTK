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

BitSeqGGMN::BitSeqGGMN(const BitStr32& bstr, size_t factor) : bstr(bstr) {
	n = bstr.length();
	ones = bstr.count();
	if(factor == 0)
		factor = bits(n - 1);
	this->factor = factor;
	b = sizeof(data_type) * Wb;
	s = b * factor;
	nRs = numSuperBlocks() + 1;
	buildRank();
}

BitSeqGGMN::BitSeqGGMN(BitStr32&& bstr, size_t factor) {
	this->bstr = bstr;
	n = bstr.length();
	ones = bstr.count();
	if(factor == 0)
		factor = bits(n - 1);
	this->factor = factor;
	b = sizeof(data_type) * Wb;
	s = b * factor;
	nRs = numSuperBlocks() + 1;
	buildRank();
}

void BitSeqGGMN::buildRank() {
	Rs.resize(nRs);
	Rs[0] = 0;
	for (size_t i = 1; i <= nRs; ++i)
		Rs[i] = Rs[i - 1] + buildRank((i-1) * factor, factor);
}

size_t BitSeqGGMN::getBytes() const {
	return BitSeq::getBytes() + bstr.getBytes() + sizeof(nRs) + sizeof(Rs) +
			sizeof(factor) + sizeof(b) + sizeof(s) + sizeof(this);
}

size_t BitSeqGGMN::rank1(size_t i) const {
	assert(bstr.getWid() == W);
	i++; /* use 1-based in rank */
	size_t r = Rs[i / s]; /* ones before this super-block */
	for(size_t b = i / s * factor; b < i / W; ++b)
		r += popcount32(bstr.getValue(b)); /* ones before this data block */
	r += popcount32(bstr.getValue(i / W) & ( (1UL << (i & lowest5_mask)) - 1)); /* ones in bits */
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
	uint32_t j = bstr.getValue(left);
	uint32_t ones = popcount32(j);
	while (ones < i) {
		i -= ones;
		left++;
		if(left > bstr.numValues())
			return length();
		j = bstr.getValue(left);
		ones = popcount32(j);
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
	uint32_t j = bstr.getValue(left);
	uint32_t zeros = W - popcount32(j);
	while(zeros < i) {
		i -= zeros;
		left++;
		if(left > bstr.numValues())
			return length();
		j = bstr.getValue(left);
		zeros = W - popcount32(j);
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
	uint32_t aux = bstr.getValue(start / W) >> des;

	if (aux) { /* has reminding */
		if(aux & 0xff)
			return start + select_tab[aux & 0xff] - 1;
		else if(aux & 0xff00)
			return start + 8 + select_tab[(aux >> 8) & 0xff] - 1;
		else if(aux & 0xff0000)
			return start + 16 + select_tab[(aux >> 16) & 0xff] - 1;
		else
			return start + 24 + select_tab[(aux >> 24) & 0xff] - 1;
	}

	/* no reminding */
	for(size_t i = start / W + 1; i < bstr.numValues(); ++i) {
		aux = bstr.getValue(i);
		if(aux) {
			if (aux & 0xff)
				return i * W + select_tab[aux & 0xff] - 1;
			else if(aux & 0xff00)
				return i * W + 8 + select_tab[(aux >> 8) & 0xff] - 1;
			else if(aux & 0xff0000)
				return i * W + 16 + select_tab[(aux >> 16) & 0xff] - 1;
			else
				return i * W + 24 + select_tab[(aux >> 24) & 0xff] - 1;
		}
	}
	return n;
}

size_t BitSeqGGMN::selectPrev1(size_t start) const {
	assert(bstr.getWid() == W);

	size_t i = start >> 5;
	size_t offset = start % W;
	uint32_t aux = bstr.getValue(i) & (- static_cast<uint32_t>(1) >> (31 - offset));

	if (aux) { /* has remaining */
		if(aux & 0xFF000000)
			return i * W + 23 + prev_tab[(aux >> 24) & 0xFF];
		else if(aux & 0xFF0000)
			return i * W + 15 + prev_tab[(aux >> 16) & 0xFF];
		else if(aux & 0xFF00)
			return i * W + 7 + prev_tab[(aux >> 8) & 0xFF];
		else
			return i * W + prev_tab[aux & 0xFF] - 1;
	}
	for(size_t k = i - 1; ; --k) {
		aux = bstr.getValue(k);
		if (aux) {
			if(aux & 0xFF000000)
				return k * W + 23 + prev_tab[(aux >> 24) & 0xFF];
			else if(aux & 0xFF0000)
				return k * W + 15 + prev_tab[(aux >> 16) & 0xFF];
			else if(aux & 0xFF00)
				return k * W + 7 + prev_tab[(aux >> 8) & 0xFF];
			else
				return k * W + prev_tab[aux & 0xFF] - 1;
		}
	}
	return 0;
}

ostream& BitSeqGGMN::save(ostream& out) const {
	BitSeq::save(out); /* save base object */
	bstr.save(out);
	out.write((const char*) &nRs, sizeof(size_t));
	out.write((const char*) Rs.c_str(), sizeof(size_t) * nRs);
	out.write((const char*) &factor, sizeof(size_t));
	out.write((const char*) &b, sizeof(size_t));
	out.write((const char*) &s, sizeof(size_t));
	return out;
}

istream& BitSeqGGMN::load(istream& in) {
	BitSeq::load(in); /* load base object */
	bstr.load(in);
	in.read((char *) &nRs, sizeof(size_t));
	size_t* tmp = new size_t[nRs];
	in.read((char*) tmp, sizeof(size_t) * nRs);
	Rs.assign(tmp, nRs);
	delete[] tmp;
	in.read((char*) &factor, sizeof(size_t));
	in.read((char*) &b, sizeof(size_t));
	in.read((char*) &s, sizeof(size_t));
	assert(nRs == numSuperBlocks() + 1);
	return in;
}

size_t BitSeqGGMN::buildRank(size_t start, size_t len) {
	size_t rank = 0;
	for(size_t i = start; i < bstr.numValues() && i < start + len; ++i)
		rank += popcount32(bstr.getValue(i));
	return rank;
}

bool operator==(const BitSeqGGMN& lhs, const BitSeqGGMN& rhs) {
	return lhs.bstr == rhs.bstr &&
			lhs.nRs == rhs.nRs && lhs.Rs == rhs.Rs &&
			lhs.factor == rhs.factor && lhs.b == rhs.b && lhs.s == rhs.s;
}

} /* namespace libSDS */
} /* namespace EGriceLab */
