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

using namespace std;

SeqIO::SeqIO(istream* in, const string& format, int maxLine) :
	in(in), out(nullptr), format(format), maxLine(maxLine) {
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
}

SeqIO::SeqIO(ostream* out, const string& format, int maxLine) :
	in(nullptr), out(out), format(format), maxLine(maxLine) {
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
}

void SeqIO::reset(istream* in, const string& format, int maxLine) {
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
	/* replace values */
	this->in = in;
	out = nullptr;
	this->format = format;
	this->maxLine = maxLine;
}

void SeqIO::reset(ostream* out, const string& format, int maxLine) {
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
	/* replace values */
	in = nullptr;
	this->out = out;
	this->format = format;
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
		throw ios_base::failure("input is not a valid FASTA format");

	*in >> name; // read the next word as id
	while(::isspace(in->peek()) && in->peek() != '\n') // ignore non-newline white spaces
		in->get();
	getline(*in, desc); // read the remaining as desc, if any
	while(in->peek() != EOF && in->peek() != fastaHead) {
		getline(*in, line);
		seq += line;
	}
	return PrimarySeq(seq, name, desc);
}

PrimarySeq SeqIO::nextFastqSeq() {
	string name, seq, desc, qual;
	char tag;
	string line;
	tag = in->get();
	if(tag != fastqHead)
		throw ios_base::failure("input is not a valid FASTQ format");
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
		for(size_t i = 0; i < seq.length(); i += maxLine)
			*out << seq.getSeq().substr(i, maxLine) << endl;
	}
	else
		*out << seq.getSeq() << endl;
}

void SeqIO::writeFastqSeq(const PrimarySeq& seq) {
	*out << fastqHead << seq.getName();
	if(!seq.getDesc().empty())
		*out << " " + seq.getDesc();
	*out << endl;
	*out << seq.getSeq() << endl;
	*out << fastqSep << endl << seq.getQual() << endl;
}

string SeqIO::guessFormat(const string& name) {
	string fn(name); /* use local copy */
	/* remove potential zip extensions */
	StringUtils::removeEnd(fn, GZIP_FILE_SUFFIX);
	StringUtils::removeEnd(fn, BZIP2_FILE_SUFFIX);
	if(StringUtils::endsWith(fn, ".fasta") || StringUtils::endsWith(fn, ".fa") ||
			StringUtils::endsWith(fn, ".fas") || StringUtils::endsWith(fn, ".fna"))
		return "fasta";
	else if(StringUtils::endsWith(fn, ".fastq") || StringUtils::endsWith(fn, ".fq"))
		return "fastq";
	else
		return "unk";
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
