/*
 * MetaGenomeAnno.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */
#include <limits>
#include <sstream>
#include "ProgEnv.h"
#include "MetaGenomeAnno.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::istringstream;

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

istream& MetaGenomeAnno::read(istream& in) {
	string line, id, name; /* current genome name */
	vector<GenomeAnno>::iterator annoIt = genomeAnnos.end();
	while(std::getline(in, line)) {
		if(StringUtils::startsWith(line, GenomeAnno::RECORD_START_TAG)) { /* a new section of genome */
			istringstream iss(line);
			GenomeAnno::readStartComment(iss, id, name);
			debugLog << "Reading GFF annotation for " << id << " (" << name << ")" << endl;
			annoIt = getAnno(id);
			if(annoIt != genomeAnnos.end()) /* a valid GenomeAnno */
				annoIt->read(in);
		}
	}
	return in;
}

ostream& MetaGenomeAnno::writeGFFHeader(ostream& out, const string& dbName, GFF::Version ver) {
	out << "##gff-version " << ver << endl;
	out << "#!processor " << progName << endl;
	out << "##metagenome " << dbName << endl;
	return out;
}

/**
 * read pre-built GFF header comments
 */
istream& MetaGenomeAnno::readGFFHeader(istream& in, string& dbName, GFF::Version& ver) {
	string tag;
	int gffVer;
	in >> tag >> gffVer;
	if(gffVer == 2)
		ver = GFF::GTF;
	else if(gffVer == 3)
		ver = GFF::GFF3;
	else
		ver = GFF::UNK;

	in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	in >> tag >> std::ws;
	std::getline(in, dbName);
	return in;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
