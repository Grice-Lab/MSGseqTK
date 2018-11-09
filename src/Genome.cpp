/*
 * Genome.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include "MSGseqTKConst.h"
#include "Genome.h"
#include "ProgEnv.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::istringstream;

const boost::regex Genome::INVALID_NAMEPREFIX_PATTERN = boost::regex("^[^\\w.:^*$@!+?-|]+");
const boost::regex Genome::INVALID_NAME_PATTERN = boost::regex("[^\\w.:^*$@!+?-|]+");
const string Genome::REPLACEMENT_STR = ".";

ostream& Genome::Chrom::save(ostream& out) const {
	StringUtils::saveString(name, out);
	out.write((const char*) &size, sizeof(uint64_t));
	return out;
}

istream& Genome::Chrom::load(istream& in) {
	StringUtils::loadString(name, in);
	in.read((char*) &size, sizeof(uint64_t));
	return in;
}

bool Genome::hasChrom(const string& chrName) const {
	for(const Chrom& chr : chroms)
		if(chr.name == chrName)
			return true;
	return false;
}

uint64_t Genome::getChromSize(const string& chrName) const {
	for(const Chrom& chr : chroms)
		if(chr.name == chrName)
			return chr.size;
	return 0;
}

size_t Genome::getChromIndex(uint64_t loc) const {
	uint64_t start = 0;
	for(vector<Chrom>::const_iterator chr = chroms.begin(); chr != chroms.end(); ++chr) {
		if(start <= loc && loc <= start + chr->size) /* include null terminal */
			return chr - chroms.begin();
		start += chr->size + 1; /* with null terminal */
	}
	return -1;
}

ostream& Genome::save(ostream& out) const {
	StringUtils::saveString(id, out);
	StringUtils::saveString(name, out);
	size_t NChrom = numChroms();
	out.write((const char*) &NChrom, sizeof(size_t));
	for(const Chrom& chr : chroms)
		chr.save(out);
	size_t Nannos = numChromAnnos();
	out.write((const char*) &Nannos, sizeof(size_t));
	for(const std::pair<string, vector<GFF>>& anno : chromAnnos) {
		StringUtils::saveString(anno.first, out);
		size_t Ngff = anno.second.size();
		out.write((const char*) &Ngff, sizeof(size_t));
		for(const GFF& gff : anno.second)
			gff.save(out);
	}
	return out;
}

istream& Genome::load(istream& in) {
	StringUtils::loadString(id, in);
	StringUtils::loadString(name, in);
	size_t NChrom = 0;
	in.read((char*) &NChrom, sizeof(size_t));
	for(size_t i = 0; i < NChrom; ++i) {
		chroms.push_back(Chrom()); /* default construct first */
		chroms[i].load(in);
	}
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

uint64_t Genome::size() const {
	uint64_t size = 0;
	for(const Chrom& chr : chroms)
		size += chr.size;
	return size + numChroms(); /* add one null gap after each chromosome */
}

bool operator==(const Genome& lhs, const Genome& rhs) {
	if(!(lhs.id == rhs.id && lhs.name == rhs.name && lhs.chroms.size() == rhs.chroms.size()))
		return false;

	for(size_t i = 0; i < lhs.chroms.size(); ++i)
		if(lhs.chroms[i] != rhs.chroms[i])
			return false;
	return true;
}

ostream& Genome::writeGFFComment(ostream& out) const {
	out << "##genome " << id << " (" << name << ")" << endl;
	return out;
}

ostream& Genome::writeGFF(ostream& out) const {
	/* write per-genome comment */
	writeGFFComment(out);
	/* write genome as first-level feature */
	GFF genomeGff(GFF_VERSION, id, progName, "genome", 1, size(), GFF::INVALID_SCORE, '.', GFF::INVALID_FRAME);
	genomeGff.setAttr("ID", id);
	genomeGff.setAttr("Name", name);
	out << genomeGff << endl;
	/* write each chromosome as second-level feature, with additional annotations */
	size_t shift = 0;
	for(const Chrom& chr : chroms) {
		UCSC::GFF chrGff(GFF_VERSION, id, progName, "chromosome", shift + 1, shift + chr.size, GFF::INVALID_SCORE, '.', GFF::INVALID_FRAME); // use genome name as seqsrc
		chrGff.setAttr("ID", chr.name);
		chrGff.setAttr("Name", chr.name);
		chrGff.setAttr("Parent", name);
		assert(chrGff.getEnd() < genomeGff.getEnd());
		out << chrGff << endl;
		/* write GFF annotations */
		if(chromAnnos.count(chr.name) > 0) { /* this chromosome has external annotations */
			for(GFF gff : chromAnnos.at(chr.name)) { /* use a local copy */
				gff.setSeqname(id); // always use genome id
				gff.shift(shift);
				if(!gff.hasAttr("Parent")) /* top level features, i.e. region, gene */
					gff.setAttr("Parent", chr.name);
				out << gff << endl;
			}
		}

		shift += chr.size + 1; /* including null terminal */
	}

	return out;
}

string Genome::formatName(const string& name) {
	return boost::replace_all_regex_copy(
			boost::replace_all_regex_copy(name, INVALID_NAMEPREFIX_PATTERN, string("")),
			INVALID_NAME_PATTERN, REPLACEMENT_STR);
}

istream& Genome::readGFF(istream& in, GFF::Version ver) {
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
