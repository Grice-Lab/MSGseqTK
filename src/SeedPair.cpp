/*
 * SMEMChain.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: zhengqi
 */

#include <algorithm>
#include "MSGseqTKConst.h"
#include "SeedPair.h"

namespace EGriceLab {
namespace MSGseqTK {

const double SeedPair::MIN_OVERRATE = 0.9;

ostream& SeedPair::save(ostream& out) const {
	out.write((const char*) &from, sizeof(int64_t));
	out.write((const char*) &to, sizeof(int64_t));
	out.write((const char*) &start, sizeof(int64_t));
	out.write((const char*) &end, sizeof(int64_t));
	out.write((const char*) &tid, sizeof(int64_t));
	out.write((const char*) &strand, sizeof(GLoc::STRAND));
	out.write((const char*) &logP, sizeof(double));
	return out;
}

istream& SeedPair::load(istream& in) {
	in.read((char*) &from, sizeof(int64_t));
	in.read((char*) &to, sizeof(int64_t));
	in.read((char*) &start, sizeof(int64_t));
	in.read((char*) &end, sizeof(int64_t));
	in.read((char*) &tid, sizeof(int64_t));
	in.read((char*) &strand, sizeof(GLoc::STRAND));
	in.read((char*) &logP, sizeof(double));
	return in;
}

ostream& SeedPair::write(ostream& out) const {
	out << from << '-' << to << ':' << start << '-' << end << ':' <<
			tid << ':' << GLoc::decode(strand) << ':' << logP;
	return out;
}

istream& SeedPair::read(istream& in) {
	in >> from;
	in.ignore(1, '-');
	in >> to;
	in.ignore(1, ':');
	in >> start;
	in.ignore(1, '-');
	in >> end;
	in.ignore(1, ':');
	in >> tid;
	in.ignore(1, ':');
	char s;
	in >> s;
	strand = GLoc::encode(s);
	in.ignore(1, ':');
	in >> logP;
	return in;
}

double SeedPair::bestLoglik(const SeedList& seeds) {
	double minLoglik = inf;
	for(const SeedPair& seed : seeds)
		minLoglik = std::min(minLoglik, seed.loglik());
	return minLoglik;
}

vector<SeedPair>& SeedPair::filter(vector<SeedPair>& pairs) {
	if(pairs.size() <= 1)
		return pairs;
	/* sort pairs by location decreasingly */
	std::sort(pairs.rbegin(), pairs.rend());
	for(vector<SeedPair>::const_iterator i = pairs.end(); i > pairs.begin(); --i) { // search backward
		for(vector<SeedPair>::const_iterator j = i - 1; j > pairs.begin(); --j) {
			if(contained(*(i - 1), *(j - 1))) { // a redundant pair
				pairs.erase(i - 1);
				break;
			}
		}
	}
	assert(!pairs.empty());
	std::reverse(pairs.begin(), pairs.end()); // reverse the order
	return pairs;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
