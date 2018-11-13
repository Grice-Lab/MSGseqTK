/*
 * GenomeAnno.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: zhengqi
 */
#include <sstream>
#include "GenomeAnno.h"

#include "ProgEnv.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::cin;
using std::cout;
using std::endl;
using std::istringstream;

const string GenomeAnno::RECORD_START_TAG = "##genome";
const string GenomeAnno::RECORD_END_TAG = "##end-genome";

size_t GenomeAnno::numAnnotated() const {
	size_t N = 0;
	for(const CHROM_ANNOMAP::value_type& entry : chromAnnos)
		N += entry.second.size();
	return N;
}

ostream& GenomeAnno::save(ostream& out) const {
	genome.save(out);
	size_t Nannos = numChromAnnotated();
	out.write((const char*) &Nannos, sizeof(size_t));
	for(const CHROM_ANNOMAP::value_type& anno : chromAnnos) {
		StringUtils::saveString(anno.first, out);
		size_t Ngff = anno.second.size();
		out.write((const char*) &Ngff, sizeof(size_t));
		for(const GFF& gff : anno.second)
			gff.save(out);
	}
	return out;
}

istream& GenomeAnno::load(istream& in) {
	genome.load(in);
	size_t Nannos = 0;
	in.read((char*) &Nannos, sizeof(size_t));
	for(size_t i = 0; i < Nannos; ++i) {
		string chr;
		StringUtils::loadString(chr, in);
		size_t Ngff = 0;
		in.read((char*) &Ngff, sizeof(size_t));
		chromAnnos[chr].resize(Ngff);
		for(size_t j = 0; j < Ngff; ++j)
			chromAnnos[chr][j].load(in);
	}
	return in;
}

ostream& GenomeAnno::write(ostream& out) const {
	const string& id = genome.getId();
	const string& name = genome.getName();
	/* write per-genome comment */
	writeStartComment(out);
	/* write genome as first-level feature */
	GFF genomeGff(FORMAT, id, progName, "genome", 1, genome.size(), GFF::INVALID_SCORE, '.', GFF::INVALID_FRAME);
	genomeGff.setAttr("ID", id);
	genomeGff.setAttr("Name", name);
	out << genomeGff << endl;

	size_t shift = 0;
	for(const Genome::Chrom& chr : genome.getChroms()) {
		GFF chrGff(FORMAT, id, progName, "chromosome", shift + 1, shift + chr.size, GFF::INVALID_SCORE, '.', GFF::INVALID_FRAME);
		chrGff.setAttr("ID", chr.name);
		chrGff.setAttr("Name", chr.name);
		assert(chrGff.getEnd() < genomeGff.getEnd());
		out << chrGff << endl;
		/* write GFF annotations */
		if(chromAnnos.count(chr.name) > 0) { /* this chromosome has external annotations */
			for(GFF gff : chromAnnos.at(chr.name)) { /* use a local copy */
				gff.setSeqname(id); // always use genome id
				gff.shift(shift);
				out << gff << endl;
			}
		}

		shift += chr.size + 1; /* including null terminal */
	}
	writeEndComment(out);
	return out;
}

ostream& GenomeAnno::writeStartComment(ostream& out) const {
	out << RECORD_START_TAG << " " << genome.getId() << " (" << genome.getName() << ")" << endl;
	return out;
}

ostream& GenomeAnno::writeEndComment(ostream& out) const {
	out << RECORD_END_TAG << endl;
	return out;
}

istream& GenomeAnno::readStartComment(istream& in, string& id, string& name) {
	string tag;
	in >> tag >> id >> name;
	name = StringUtils::stripQuotes(name, "()");
	return in;
}

istream& GenomeAnno::read(istream& in, GFF::Version ver) {
	string line;
	GFF gffRecord(ver);
	while(std::getline(in, line)) {
		if(line.empty())
			continue;
		else if(line.front() == GFF::COMMENT_CHAR) {
			if(StringUtils::startsWith(line, "##gff-version 3")) {
				gffRecord.setVer(GFF::GFF3);
				debugLog << "  GFF version determined by embedded comment" << endl;
			}
			else if(StringUtils::startsWith(line, "##gff-version 2")) {
				gffRecord.setVer(GFF::GTF);
				debugLog << "  GFF version determined by embedded comment" << endl;
			}
			else if(line == RECORD_END_TAG)
				break;
			else
				continue;
		}
		else {
			istringstream iss(line);
			iss >> gffRecord;
			const string& chr = gffRecord.getSeqname();
			chromAnnos[chr].push_back(gffRecord);
		}
	}

	return in;
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
