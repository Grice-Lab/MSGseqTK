/*
 * Genome.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#include <boost/algorithm/string.hpp>
#include <unordered_set>
#include <algorithm>
#include "MSGseqTKConst.h"
#include "Genome.h"
#include "ProgEnv.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::unordered_set;
const std::regex Genome::INVALID_NAMEPREFIX_PATTERN = std::regex("^[^\\w.:^*$@!+?-|]+");
const std::regex Genome::INVALID_NAME_PATTERN = std::regex("[^\\w.:^*$@!+?-|]+");
const string Genome::REPLACEMENT_STR = ".";

ostream& Genome::Chrom::save(ostream& out) const {
	StringUtils::saveString(name, out);
	out.write((const char*) &len, sizeof(int64_t));
	return out;
}

ostream& Genome::Chrom::saveSeq(ostream& out) const {
	StringUtils::saveString(seq, out, len + 1); // include null terminal
	return out;
}

istream& Genome::Chrom::load(istream& in) {
	StringUtils::loadString(name, in);
	in.read((char*) &len, sizeof(int64_t));
	return in;
}

istream& Genome::Chrom::loadSeq(istream& in) {
	StringUtils::loadString(seq, in, len);
	assert(0 == in.get());
	return in;
}

size_t Genome::removeRedundantChroms() {
	size_t n = numChroms(); /* old number of chromosomes */
	unordered_set<string> names; /* unique chrom name set */
	/* remove-erase redundnt chromsomes */
	chroms.erase(std::remove_if(chroms.begin(), chroms.end(),
			[&] (const Chrom& chr)->bool { return ! names.insert(chr.name).second; }),
			chroms.end());
	return n - numChroms();
}

size_t Genome::getChromIndex(int64_t loc) const {
	int64_t start = 0;
	for(vector<Chrom>::const_iterator chr = chroms.begin(); chr != chroms.end(); ++chr) {
		if(start <= loc && loc <= start + chr->size()) /* include null terminal */
			return chr - chroms.begin();
		start += chr->size() + 1; /* with null terminal */
	}
	return -1;
}

ostream& Genome::save(ostream& out) const {
	StringUtils::saveString(id, out);
	StringUtils::saveString(name, out);
	size_t NC = numChroms();
	out.write((const char*) &NC, sizeof(size_t));
	for(Chrom chr : chroms)
		chr.save(out);
	return out;
}

ostream& Genome::saveSeq(ostream& out) const {
	for(const Chrom& chr : chroms)
		chr.saveSeq(out);
	return out;
}

istream& Genome::load(istream& in) {
	StringUtils::loadString(id, in);
	StringUtils::loadString(name, in);
	size_t NC = 0;
	in.read((char*) &NC, sizeof(size_t));
	chroms.resize(NC);
	for(size_t i = 0; i < NC; ++i)
		chroms[i].load(in);
	return in;
}

istream& Genome::loadSeq(istream& in) {
	for(Chrom& chr : chroms)
		chr.loadSeq(in);
	return in;
}

int64_t Genome::size() const {
	int64_t len = 0;
	for(const Chrom& chr : chroms)
		len += chr.size();
	return len;
}

string Genome::formatName(const string& name) {
	return std::regex_replace(
			std::regex_replace(name, INVALID_NAMEPREFIX_PATTERN, ""),
			INVALID_NAME_PATTERN, REPLACEMENT_STR);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
