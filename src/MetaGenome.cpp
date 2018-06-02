/*
 * Genome.cpp
 *
 *  Created on: May 4, 2018
 *      Author: zhengqi
 */

#include "MetaGenome.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqClean {
using Eigen::Map;

MetaGenome& MetaGenome::addGenome(const string& genomeName, const DNAseq& genomeSeq) {
	genomeNames.push_back(genomeName);
	size += genomeSeq.length();
	for(DNAseq::value_type b : genomeSeq)
		if(DNAalphabet::isBase(b))
			baseCount[b - DNAalphabet::A]++;

	return *this;
}

ostream& MetaGenome::save(ostream& out) const {
	/* save basic info */
	out.write((const char*) &size, sizeof(uint64_t));
	/* save chromosomes */
	size_t NGenome = genomeNames.size();
	out.write((const char*) &NGenome, sizeof(size_t));
	for(const string& chr : genomeNames)
		StringUtils::saveString(chr, out);
	/* save base counts */
	uint64_t bufi[4];
	Map<Vector4l> baseCountMap(bufi);
	baseCountMap = baseCount; /* copy data */
	out.write((const char*) bufi, sizeof(uint64_t) * 4);
	return out;
}

istream& MetaGenome::load(istream& in) {
	/* load basic info */
	in.read((char*) &size, sizeof(uint64_t));
	/* load chromosomes */
	size_t NGenome = 0;
	in.read((char *) &NGenome, sizeof(size_t));
	for(size_t i = 0; i < NGenome; ++i)
		StringUtils::loadString(genomeNames[i], in);
	/* load base counts */
	uint64_t bufi[4];
	in.read((char*) bufi, sizeof(uint64_t) * 4);
	baseCount = Map<Vector4l>(bufi); /* copy data */
	return in;
}

MetaGenome& MetaGenome::operator+=(const MetaGenome& other) {
	genomeNames.insert(genomeNames.end(), other.genomeNames.begin(), other.genomeNames.end());
	size += other.getSize();
	baseCount += other.baseCount;
	return *this;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
