/*
 * MetaGenomeAnno.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */

#include "MetaGenomeAnno.h"

namespace EGriceLab {
namespace MSGseqTK {

ostream& MetaGenomeAnno::save(ostream& out) const {
	size_t N = numAnnotated();
	out.write((const char*) &N, sizeof(size_t));
	for(const vector<GenomeAnno>::value_type& anno : genomeAnnos)
		anno.save(out);
	return out;
}

istream& MetaGenomeAnno::load(istream& in) {
	size_t N = 0;
	in.read((char*) &N, sizeof(size_t));
	genomeAnnos.resize(N);
	for(size_t i = 0; i < N; ++i)
		genomeAnnos[i].load(in);
	return in;
}

ostream& MetaGenomeAnno::write(ostream& out) const {
	/* write each genome with external annotations */
	for(const GenomeAnno& genome : genomeAnnos)
		genome.write(out);
	return out;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
