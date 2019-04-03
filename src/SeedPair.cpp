/*
 * SMEMChain.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: zhengqi
 */

#include "MSGseqTKConst.h"
#include "SeedPair.h"

namespace EGriceLab {
namespace MSGseqTK {

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

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
