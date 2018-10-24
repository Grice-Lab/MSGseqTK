/*
 * Seq.cpp
 *
 *  Created on: Oct 17, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "Seq.h"

namespace EGriceLab {
namespace libSDS {

size_t Seq::access(size_t i, size_t & r) const {
	size_t s = access(i);
	r = rank(s, i);
	return s;
}

size_t Seq::rank(size_t s, size_t i) const {
	size_t r = 0;
	for(size_t k = 0; k <= i && k < n; ++k)
		if(access(k) == s)
			r++;
	return r;
}

size_t Seq::select(size_t s, size_t r) const {
	for(size_t k = 0, count = 0; k < n; ++k) {
		if(access(k) == s) {
			count++;
			if(count == r)
				return k;
		}
	}
	return n;
}

size_t Seq::selectNext(size_t s, size_t start) const {
	size_t r = rank(s, start);
	return select(s, r + 1);
}

size_t Seq::maxValue() const {
	size_t max = 0;
	for(size_t k = 0; k < n; ++k)
		max = std::max(max, access(k));
	return max;
}

ostream& Seq::save(ostream& out) const {
	out.write((const char*) &n, sizeof(size_t));
	out.write((const char*) &sigma, sizeof(size_t));
	return out;
}

istream& Seq::load(istream& in) {
	in.read((char*) &n, sizeof(size_t));
	in.read((char*) &sigma, sizeof(size_t));
	return in;
}

} /* namespace libSDS */
} /* namespace EGriceLab */
