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
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

using namespace std;
const string MetaGenomeAnno::GENOME_START_TAG = "##genome";
const string MetaGenomeAnno::GENOME_END_TAG = "##end-genome";
const string MetaGenomeAnno::METAGENOME_ID_TAG = "metagenomeID";
const string MetaGenomeAnno::METAGENOME_NAME_TAG = "metagenomeName";
const string MetaGenomeAnno::METAGENOME_TYPE_TAG = "metagenomeType";
const string MetaGenomeAnno::EXTERNAL_ID_TAG = "ID";
const string MetaGenomeAnno::EXTERNAL_NAME_TAG = "Name";


size_t MetaGenomeAnno::writeGenomeAnnos(ostream& out, const Genome& genome, const vector<GFF>& gffRecords) {
	size_t n = 0;
	/* start per-genome comment */
	writeStartComment(out, genome);

	/* write genome annotation */
	out << getAnno(genome) << endl;
	n++;

	/* write chrom annotations */
	for(const Genome::Chrom& chr : genome.chroms)
		out << getAnno(genome, chr) << endl;
	n += genome.numChroms();

	/* write auxilary annotations */
	for(const GFF& gff : gffRecords)
		out << gff << endl;
	n += gffRecords.size();

	/* end comment */
	writeEndComment(out);
	return n;
}

vector<GFF> MetaGenomeAnno::read(istream& in, GFF::Version ver) {
	string line;
	vector<GFF> gffRecords;
	while(std::getline(in, line)) {
		if(line.empty())
			continue;
		else if(line.front() == GFF::COMMENT_CHAR) {
			if(StringUtils::startsWith(line, "##gff-version 3")) {
				ver = GFF::GFF3;
				debugLog << "  GFF version determined by embedded comment" << endl;
			}
			else if(StringUtils::startsWith(line, "##gff-version 2")) {
				ver = GFF::GTF;
				debugLog << "  GFF version determined by embedded comment" << endl;
			}
			else if(line == "##FASTA") /* appended FASTA sequences */
				break;
			else
				continue;
		}
		else {
			// read actual records
			if(ver == GFF::UNK) {
				cerr << "  Unable to guess GFF version from file content" << endl;
				in.setstate(std::ios_base::badbit);
				break;;
			}
			istringstream iss(line);
			GFF gff;
			if(gff.read(iss, ver).bad())
				continue;
			else
				gffRecords.push_back(gff);
		}
	}
	return gffRecords;
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

istream& MetaGenomeAnno::readStartComment(istream& in, Genome& genome) {
	string tag;
	in >> tag >> genome.id >> genome.name;
	genome.name = StringUtils::stripQuotes(genome.name, "()");
	return in;
}

ostream& MetaGenomeAnno::writeStartComment(ostream& out, const Genome& genome) {
	out << GENOME_START_TAG << " " << genome.displayId() << endl;
	return out;
}

ostream& MetaGenomeAnno::writeEndComment(ostream& out) {
	out << GENOME_END_TAG << endl;
	return out;
}

vector<GFF>& MetaGenomeAnno::addMetagenomeTags(const Genome& genome, vector<GFF>& gffRecords,
		const string& idTag, const string& nameTag) {
	string id;
	string name;
	string type;
	for(GFF& record : gffRecords) {
		if(!record.hasAttr("Parent")) { //update only when it is not a child feature
			id = genome.getId() + ":" + record.getAttr(idTag);
			name = record.hasAttr(nameTag) ? genome.getName() + ":" + record.getAttr(nameTag) : id;
			type = record.getType();
		}
		record.setAttr(METAGENOME_ID_TAG, id);
		record.setAttr(METAGENOME_NAME_TAG, name);
		record.setAttr(METAGENOME_TYPE_TAG, type);
	}
	return gffRecords;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
