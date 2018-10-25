/*
 * SeqRRR.cpp
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "WaveletTreeRRR.h"

namespace EGriceLab {
namespace libSDS {

ostream& WaveletTreeRRR::save(ostream& out) const {
	/* save base class fields */
	Seq::save(out);
	out.write((const char*) &height, sizeof(size_t));
	out.write((const char*) &wid, sizeof(size_t));
	out.write((const char*) &min, sizeof(size_t));
	out.write((const char*) &max, sizeof(size_t));
	size_t nOCC = OCC.size();
	out.write((const char*) &nOCC, sizeof(size_t));
	out.write((const char*) OCC.data(), sizeof(size_t) * nOCC);
	/* write all bseqs */
	for(const BitSeqRRR& bs : bseqs)
		bs.save(out);
	return out;
}

istream& WaveletTreeRRR::load(istream& in) {
	/* load base class fields */
	Seq::load(in);
	in.read((char*) &height, sizeof(size_t));
	in.read((char*) &wid, sizeof(size_t));
	in.read((char*) &min, sizeof(size_t));
	in.read((char*) &max, sizeof(size_t));
	size_t nOCC;
	in.read((char*) &nOCC, sizeof(size_t));
	OCC.resize(nOCC);
	in.read((char*) OCC.data(), sizeof(size_t) * nOCC);
	/* read all baseqs */
	bseqs.resize(height);
	for(size_t i = 0; i < height; ++i)
		bseqs[i].load(in);
	return in;
}

size_t WaveletTreeRRR::access(size_t i) const {
	size_t s = 0;
	size_t start = 0;
	for (size_t level = 0; level < height; ++level) {
		size_t optR = 0;
		size_t before = 0;
		if (start > 0)
			before = bseqs[level].rank1(start - 1);
		if (bseqs[level].access(i, optR)) {
			s |= 1UL << (height - level - 1);
			i = optR - 1 - before;
			start = OCC[s];
			i += start;
		}
		else
			i = optR - 1 + before;
	}
	return s;
}

size_t WaveletTreeRRR::access(size_t i, size_t& r) const {
	size_t s = 0;
	size_t start = 0;
	for (size_t level = 0; level < height; ++level) {
		size_t optR = 0;
		size_t before = 0;
		if(start > 0)
			before = bseqs[level].rank1(start - 1);

		if(bseqs[level].access(i, optR)) {
			s |= (1UL << (height - level - 1));
			r = optR - before;
			start = OCC[s];
			i = r - 1 + start;
		}
		else {
			i = optR - 1 + before;
			r = i + 1 - start;
		}
	}
	return s;
}

size_t WaveletTreeRRR::rank(size_t s, size_t i) const {
	if(i >= n)
		return -1;
	size_t start = 0;
	size_t r = 0;

	for(size_t level = 0; level < height; ++level) {
		size_t masked = (s >> (height - level - 1)) << (height - level - 1);
		size_t before = 0;
		if (start > 0)
			before = bseqs[level].rank1(start - 1);

		if (test(s, level)) {
			r = bseqs[level].rank1(i) - before;
			start = OCC[masked];
			i = r + start - 1;
		}
		else {
			r = i - start + before - bseqs[level].rank1(i) + 1;
			masked += (1 << (height - level - 1));
			i = r + start - 1;
		}
		if(r == 0)
			break;
	}
	return r;
}

size_t WaveletTreeRRR::select(size_t s, size_t r) const {
	size_t mask = (1UL << height) - 2;
	size_t sum = 2;
	size_t i = r;

	for(size_t level = height; level > 0; --level) {
		size_t start = get_start(s, mask);
		start = OCC[start];

		size_t ones_start = 0;
		if (start > 0)
			ones_start = bseqs[level - 1].rank1(start - 1);

		if (test(s, level - 1))
			i = bseqs[level - 1].select1(ones_start + i) - start + 1;
		else
			i = bseqs[level - 1].select0(start - ones_start + i) - start + 1;

		mask <<= 1;
		sum <<= 1;
	}
	return i - 1;
}

} /* namespace libSDS */
} /* namespace EGriceLab */
