/*
 * SeqIO.h
 *
 *  Created on: Apr 25, 2018
 *      Author: zhengqi
 */

#ifndef SRC_SEQIO_H_
#define SRC_SEQIO_H_
#include <string>
#include <iostream>
#include <stdexcept>
#include "PrimarySeq.h"

namespace EGriceLab {
namespace MSGseqTK {

using std::string;
using std::istream;
using std::ostream;
using std::streambuf;
using std::ifstream;
using std::ofstream;

/**
 * A class to handle IO operation for PrimarySeq of various format
 */
class SeqIO {
public:
	/* nested types and enums */
	enum FORMAT { UNK, FASTA, FASTQ };

	/* constructors */
	/** disable default constructor */
	SeqIO() = default;

	/* Disable copy and assign constructors */
	SeqIO(const SeqIO& other) = delete;
	SeqIO& operator=(const SeqIO& other) = delete;

	/**
	 * Construct a SeqIO object in READ mode with given info
	 */
	SeqIO(istream* in, FORMAT fmt, bool fixQual = false, int maxLine = DEFAULT_MAX_LINE);

	/**
	 * Construct a SeqIO object in WRITE mode with given info
	 */
	SeqIO(ostream* out, FORMAT fmt, bool fixQual = false, int maxLine = DEFAULT_MAX_LINE);

public:
	/* Getters and Setters */
	FORMAT getFormat() const {
		return fmt;
	}

	int getMaxLine() const {
		return maxLine;
	}

	void setMaxLine(int maxLine) {
		this->maxLine = maxLine;
	}

	/* member methods */
	/** set the input to a given a new istream, will not close the old one */
	void reset(istream* in, FORMAT fmt, bool fixQual = false, int maxLine = DEFAULT_MAX_LINE);

	/** set the out to a given a new ostream, will not close the old one */
	void reset(ostream* out, FORMAT fmt, bool fixQual = false, int maxLine = DEFAULT_MAX_LINE);

	/**
	 * test whether this file has next PrimarySeq
	 * @return true if everything is good and has symbol indicating nextSeq exists
	 */
	bool hasNext();

	/**
	 * Get next PrimarySeq, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextSeq not available or other IO exception
	 * @throw
	 */
	PrimarySeq nextSeq();

	/**
	 * Write a seq to the output
	 * @param seq  a PrimarySeq
	 * @throw std::ios_base::failure if any IO exception
	 */
	void writeSeq(const PrimarySeq& seq);

private:
	/**
	 * Get next PrimarySeq in fasta format, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextSeq not available or other IO exception
	 */
	PrimarySeq nextFastaSeq();

	/**
	 * Get next PrimarySeq in fasta format, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextSeq not available or other IO exception
	 */
	PrimarySeq nextFastqSeq();

	/**
	 * test whether this file has next PrimarySeq in fasta format
	 * @return true if everything is good and has symbol indicating nextSeq exists
	 */
	bool hasNextFasta();

	/**
	 * test whether this file has next PrimarySeq in fastq format
	 * @return true if everything is good and has symbol indicating nextSeq exists
	 */
	bool hasNextFastq();

	/**
	 * Write a seq to the output in fasta format
	 * @param seq  a PrimarySeq
	 * @throw std::ios_base::failure if any IO exception
	 */
	void writeFastaSeq(const PrimarySeq& seq);

	/**
	 * Write a seq to the output in fastq format,
	 * with maxLine restricted
	 * @param seq  a PrimarySeq
	 * @param maxLine  max characters in a line, set to -1 for limits
	 * @throw std::ios_base::failure if any IO exception
	 */
	void writeFastqSeq(const PrimarySeq& seq);

private:
	/** member fields */
	FORMAT fmt;
	int maxLine;
	bool fixQual = false;

	istream* in; /* input */
	ostream* out; /* output */

	/* static members */
public:
	static const char fastaHead = '>';
	static const char fastqHead = '@';
	static const int DEFAULT_MAX_LINE = 60;
	static const char fastqSep = '+';

	/** guess sequence file format by its extension */
	static FORMAT guessFormat(const string& name);
};

inline bool SeqIO::hasNext() {
	switch(fmt) {
	case FASTA:
		return hasNextFasta();
	case FASTQ:
		return hasNextFastq();
	default:
		return false;
	}
}

inline PrimarySeq SeqIO::nextSeq() {
	switch(fmt) {
	case FASTA:
		return nextFastaSeq();
	case FASTQ:
		return nextFastqSeq();
	default:
		throw std::ios_base::failure("Unsupported sequence format");
	}
}

inline void SeqIO::writeSeq(const PrimarySeq& seq) {
	switch(fmt) {
	case FASTA:
		return writeFastaSeq(seq);
	case FASTQ:
		return writeFastqSeq(seq);
	default:
		throw std::ios_base::failure("Unsupported sequence format");
	}
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_SEQIO_H_ */
