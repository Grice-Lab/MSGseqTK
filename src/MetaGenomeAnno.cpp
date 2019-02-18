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

size_t MetaGenomeAnno::numAnnotations() const {
	size_t N = 0;
	for(const GENOME_ANNOMAP::value_type& entry : genomeAnnos)
		N += entry.second.size();
	for(const CHROM_ANNOMAP::value_type& entry : chromAnnos)
		N += entry.second.size();
	return N;
}

void MetaGenomeAnno::addGenome(const Genome& genome) {
	/* add a genome-level annotation */
	addGenomeAnno(genome.id, GFF(genome.id, progName, "genome",
			1, genome.size(), GFF::INVALID_SCORE, '.', GFF::INVALID_FRAME,
			GFF::attr_map { {"ID", genome.id}, {"Name", Genome::formatName(genome.name)} }));

	for(const Genome::Chrom& chr : genome.chroms) {
		string chrId = MetaGenome::getChromId(genome.id, chr.name);
		addChromAnno(chrId, GFF(chr.name, progName, "chromosome",
				1, chr.size(), GFF::INVALID_SCORE, '.', GFF::INVALID_FRAME,
				GFF::attr_map { {"ID", chrId}, {"Name", chr.name}, {"Parent", genome.id} }));
	}
}

ostream& MetaGenomeAnno::save(ostream& out) const {
	/* write annotations */
	const size_t NG = numAnnotatedGenomes();
	const size_t NC = numAnnotatedChroms();
	out.write((const char*) &NG, sizeof(size_t));
	out.write((const char*) &NC, sizeof(size_t));

	for(const GENOME_ANNOMAP::value_type& entry : genomeAnnos) {
		StringUtils::saveString(entry.first, out);
		const size_t NA = entry.second.size();
		out.write((const char*) &NA, sizeof(size_t));
		for(const GFF& gff : entry.second)
			gff.save(out);
	}

	for(const CHROM_ANNOMAP::value_type& entry : chromAnnos) {
		StringUtils::saveString(entry.first, out);
		const size_t NA = entry.second.size();
		out.write((const char*) &NA, sizeof(size_t));
		for(const GFF& gff : entry.second)
			gff.save(out);
	}

	return out;
}

istream& MetaGenomeAnno::load(istream& in) {
	/* load annotations */
	size_t NG = 0;
	size_t NC = 0;
	in.read((char*) &NG, sizeof(size_t));
	in.read((char*) &NC, sizeof(size_t));

	for(size_t i = 0; i < NG; ++i) {
		GENOME_ANNOMAP::key_type key;
		StringUtils::loadString(key, in);
		size_t NA = 0;
		in.read((char*) &NA, sizeof(size_t));
		genomeAnnos[key].resize(NA); // default construct
		for(size_t i = 0; i < NA; ++i)
			genomeAnnos[key][i].load(in);
	}

	for(size_t i = 0; i < NC; ++i) {
		CHROM_ANNOMAP::key_type key;
		StringUtils::loadString(key, in);
		size_t NA = 0;
		in.read((char*) &NA, sizeof(size_t));
		chromAnnos[key].resize(NA); // default construct
		for(size_t i = 0; i < NA; ++i)
			chromAnnos[key][i].load(in);
	}

	return in;
}

ostream& MetaGenomeAnno::write(ostream& out, const Genome& genome) const {
	/* start per-genome comment */
	writeStartComment(out, genome);

	/* write genome-level annotations, if any */
	if(genomeAnnos.count(genome.id) > 0) {
		for(const GFF& genomeGff : genomeAnnos.at(genome.id))
			out << genomeGff << endl;
	}

	/* write chrom-level annotations, if any */
	for(const Genome::Chrom& chr : genome.chroms) {
		string chrId = MetaGenome::getChromId(genome.id, chr.name);
		if(chromAnnos.count(chrId) > 0) {
			for(const GFF& chrGff : chromAnnos.at(chrId))
				out << chrGff << endl;
		}
	}

	/* end per-genome comment */
	writeEndComment(out);
	return out;
}

istream& MetaGenomeAnno::read(istream& in, const Genome& genome, GFF::Version ver) {
	string line;
	GFF gffRecord;
	Genome gffGenome;
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
			else if(StringUtils::startsWith(line, GENOME_START_TAG)) {
				readStartComment(in, gffGenome);
				if(!(gffGenome.id == genome.id && gffGenome.name == genome.name)) {
					warningLog << "  GFF file content doesn't match given genome information" << endl;
					in.setstate(std::ios_base::badbit);
					break;
				}
			}
			else if(StringUtils::startsWith(line, GENOME_END_TAG))
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
			if(gffRecord.read(iss, ver).bad())
				break;
			if(gffRecord.getSource() == progName) // genome-level anno
				addGenomeAnno(genome.id, gffRecord);
			else
				addChromAnno(MetaGenome::getChromId(genome.id, gffRecord.getSeqname()), gffRecord);
		}
	}

	return in;
}

ostream& MetaGenomeAnno::write(ostream& out, const MetaGenome& mtg) const {
	/* write genome annotations */
	for(const Genome& genome : mtg.getGenomes())
		write(out, genome);
	return out;
}

istream& MetaGenomeAnno::read(istream& in, const MetaGenome& mtg) {
	string line;
	Genome genome;
	while(std::getline(in, line)) {
		if(line.empty())
			continue;
		else if(line.front() == GFF::COMMENT_CHAR) {
			if(StringUtils::startsWith(line, GENOME_START_TAG)) { /* a new section of genome */
				istringstream iss(line);
				readStartComment(iss, genome);
				if(!mtg.hasGenome(genome.id)) {
					cerr << "Error: GFF file contains genome " << genome.displayId() << " that is not found in the MetaGenome" << endl;
					in.setstate(std::ios_base::badbit);
					break;
				}
				debugLog << "Reading GFF annotation for " << genome.displayId() << endl;
				read(in, genome); // read till next genome
			}
		}
		else
			continue;
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

MetaGenomeAnno& MetaGenomeAnno::operator+=(const MetaGenomeAnno& other) {
	return *this;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
