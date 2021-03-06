/*
 * SeqIO.cpp
 *
 *  Created on: Apr 25, 2018
 *      Author: zhengqi
 */

#include <fstream>
#include "SeqIO.h"
#include "MSGseqTKConst.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace MSGseqTK {
using std::endl;

void SeqIO::reset(istream* in, FORMAT fmt, bool maskLower, bool fixQual, int maxLine) {
	/* check format support */
	if(fmt == UNK)
		throw std::invalid_argument("Unsupported file format");
	/* replace values */
	this->in = in;
	out = nullptr;
	this->fmt = fmt;
	this->maskLower = maskLower;
	this->fixQual = fixQual;
	this->maxLine = maxLine;
}

void SeqIO::reset(ostream* out, FORMAT fmt, bool maskLower, bool fixQual, int maxLine) {
	/* check format support */
	if(fmt == UNK)
		throw std::invalid_argument("Unsupported file format");
	/* replace values */
	in = nullptr;
	this->out = out;
	this->fmt = fmt;
	this->maskLower = maskLower;
	this->fixQual = fixQual;
	this->maxLine = maxLine;
}

bool SeqIO::hasNextFasta() {
	char c = in->peek();
	return c != EOF && c == fastaHead;
}

bool SeqIO::hasNextFastq() {
	char c = in->peek();
	return c != EOF && c == fastqHead;
}

PrimarySeq SeqIO::nextFastaSeq() {
	string name, seq, desc;
	char tag;
	string line;
	tag = in->get();
	if(tag != fastaHead)
		throw std::ios_base::failure("input is not a valid FASTA format");

	*in >> name; // read the next word as id
	while(::isspace(in->peek()) && in->peek() != '\n') // ignore non-newline white spaces
		in->get();
	getline(*in, desc); // read the remaining as desc, if any
	while(in->peek() != EOF && in->peek() != fastaHead) {
		getline(*in, line);
		seq += line;
	}
	if(maskLower)
		dna::mask(seq);
	return !fixQual ? PrimarySeq(seq, name, desc) : PrimarySeq(seq, name, desc).fixQual();
}

PrimarySeq SeqIO::nextFastqSeq() {
	string name, seq, desc, qual;
	char tag;
	string line;
	tag = in->get();
	if(tag != fastqHead)
		throw std::ios_base::failure("input is not a valid FASTQ format");
	*in >> name; // read the next word as id
	while(::isspace(in->peek()) && in->peek() != '\n') // ignore non-newline white spaces
		in->get();
	getline(*in, desc); // read the remaining as description
	getline(*in, seq);  // read seq line
	getline(*in, line); // ignore sep line
	getline(*in, qual); // read qual line
	return PrimarySeq(seq, name, desc, qual);
}

void SeqIO::writeFastaSeq(const PrimarySeq& seq) {
	*out << fastaHead << seq.getName();
	if(!seq.getDesc().empty())
		*out << " " << seq.getDesc();
	*out << endl;
	if(maxLine > 0) {
		for(size_t i = 0; i < seq.length(); i += maxLine) {
			dna::write(seq.getSeq().substr(i, maxLine), *out);
			*out << endl;
		}
	}
	else {
		dna::write(seq.getSeq(), *out);
		*out << endl;
	}
}

void SeqIO::writeFastqSeq(const PrimarySeq& seq) {
	*out << fastqHead << seq.getName();
	if(!seq.getDesc().empty())
		*out << " " + seq.getDesc();
	*out << endl;
	dna::write(seq.getSeq(), *out);
	*out << endl << fastqSep << endl;
	quality::write(seq.getQual(), *out);
	*out << endl;
}

SeqIO::FORMAT SeqIO::guessFormat(const string& name) {
	string fn(name); /* use local copy */
	/* remove potential zip extensions */
	StringUtils::removeEnd(fn, GZIP_FILE_SUFFIX);
	StringUtils::removeEnd(fn, BZIP2_FILE_SUFFIX);
	if(StringUtils::endsWith(fn, ".fasta") || StringUtils::endsWith(fn, ".fa") ||
			StringUtils::endsWith(fn, ".fas") || StringUtils::endsWith(fn, ".fna"))
		return FASTA;
	else if(StringUtils::endsWith(fn, ".fastq") || StringUtils::endsWith(fn, ".fq"))
		return FASTQ;
	else
		return UNK;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
