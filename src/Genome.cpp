/*
 * Genome.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#include "Genome.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqClean {

using namespace Eigen;

Genome::Genome(const string& name, const DNAseq& seq)
: name(name) {
	size = seq.length();
	for(DNAseq::value_type b : seq)
		if(DNAalphabet::isBase(b))
			baseCount[b - DNAalphabet::A]++;
}

ostream& Genome::save(ostream& out) const {
	StringUtils::saveString(name, out);
	out.write((const char*) &size, sizeof(uint64_t));
	uint64_t bufi[4];
	Map<Vector4l> baseCountMap(bufi);
	baseCountMap = baseCount; /* copy data */
	out.write((const char*) bufi, 4 * sizeof(uint64_t));
	return out;
}

istream& Genome::load(istream& in) {
	StringUtils::loadString(name, in);
	in.read((char*) &size, sizeof(uint64_t));
	uint64_t bufi[4];
	in.read((char*) bufi, 4 * sizeof(uint64_t));
	baseCount = Map<Vector4l>(bufi); /* copy data */
	return in;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
