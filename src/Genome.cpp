/*
 * Genome.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: zhengqi
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include "MSGseqTKConst.h"
#include "Genome.h"
#include "ProgEnv.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {

const boost::regex Genome::INVALID_NAMEPREFIX_PATTERN = boost::regex("^[^\\w.:^*$@!+?-|]+");
const boost::regex Genome::INVALID_NAME_PATTERN = boost::regex("[^\\w.:^*$@!+?-|]+");
const string Genome::REPLACEMENT_STR = ".";

ostream& Genome::Chrom::save(ostream& out) const {
	StringUtils::saveString(name, out);
	seq.nt16Save(out); /* use compressed saving */
	return out;
}

istream& Genome::Chrom::load(istream& in) {
	StringUtils::loadString(name, in);
	seq.nt16Load(in);
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
			return chr.size();
	return 0;
}

ostream& Genome::save(ostream& out) const {
	StringUtils::saveString(id, out);
	StringUtils::saveString(name, out);
	size_t NC = numChroms();
	out.write((const char*) &NC, sizeof(size_t));
	for(const Chrom& chr : chroms)
		chr.save(out);
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

uint64_t Genome::size() const {
	uint64_t size = 0;
	for(const Chrom& chr : chroms)
		size += chr.size();
	return size + numChroms(); /* add one null gap after each chromosome */
}

DNAseq Genome::getSeq() const {
	DNAseq seq;
	seq.reserve(size() - 1);
	for(const Chrom& chr : chroms) {
		if(!seq.empty())
			seq.push_back(DNAalphabet::GAP_BASE);
		seq += chr.seq;
	}
	return seq;
}

DNAseq Genome::getSeqRevOrder() const {
	cerr << "getRevOrder for: " << id << endl;
	DNAseq revSeq;
	revSeq.reserve(size() - 1);
	for(vector<Chrom>::const_reverse_iterator chr = chroms.rbegin(); chr != chroms.rend(); ++chr) {
		cerr << "  addingChrom RevOrder for: " << chr->name << endl;
		if(!revSeq.empty())
			revSeq.push_back(DNAalphabet::GAP_BASE);
		revSeq += chr->seq;
	}
	return revSeq;
}

string Genome::formatName(const string& name) {
	return boost::replace_all_regex_copy(
			boost::replace_all_regex_copy(name, INVALID_NAMEPREFIX_PATTERN, string("")),
			INVALID_NAME_PATTERN, REPLACEMENT_STR);
}

} /* namespace MSGseqTK */
} /* namespace EGriceLab */
