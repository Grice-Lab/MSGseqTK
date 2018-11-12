/*
 * MetaGenomeAnno.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */

#include "MetaGenomeAnno.h"

namespace EGriceLab {
namespace MSGseqTK {

ostream& MetaGenomeAnno::write(ostream& out) const {
	/* write each genome with external annotations */
	for(const GenomeAnno& genome : genomeAnnos)
		genome.write(out);
	return out;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
