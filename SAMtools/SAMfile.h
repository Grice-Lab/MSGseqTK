/*
 * SAMfile.h
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */

#ifndef SRC_SAMFILE_H_
#define SRC_SAMFILE_H_

#include <string>
#include <cstdint>
#include <utility>
#include <htslib/sam.h>

#include "BAM.h"
#include "BAMheader.h"
#include "HTSindex.h"

namespace EGriceLab {
namespace SAMtools {

using std::string;

/**
 * A C++ wrapper class for the latest SAM file,
 * derived from with the samFile/htsFile struct described from htslib C library
 */
class SAMfile {
public:
	/* constructors */
	/** construct a SAMfile from a given path and mode */
	SAMfile(const string& path, const string& mode) {
		samFh = sam_open(path.c_str(), mode.c_str());
	}

	/** no copy constructor */
	SAMfile(const SAMfile&) = delete;
	/** no copy assignment */
	SAMfile& operator=(const SAMfile&) = delete;

	/** move constructor */
	SAMfile(SAMfile&&) = default;

	/** move assignment */
	SAMfile& operator=(SAMfile&& other) {
		sam_close(samFh);
		samFh = std::move(other.samFh);
		header = std::move(other.header);
		idx = std::move(other.idx);
		return *this;
	}

	/** destructor */
	virtual ~SAMfile() {
		sam_close(samFh);
	}

	/* member methods */
	/** get BAMheader of this file */
	const BAMheader& getHeader() const {
		return header;
	}

	/** set header of this file */
	void setHeader(const BAMheader& header) {
		this->header = header;
	}

	/**
	 * load an existing index from a given BAM file basename and optional suffix
	 */
	SAMfile& loadIndex(const string& basename, const string& suffix = HTSindex::DEFAULT_INDEX_SUFFIX) {
		if(suffix.empty())
			idx = HTSindex(sam_index_load(samFh, basename.c_str()));
		else
			idx = HTSindex(sam_index_load2(samFh, basename.c_str(), (basename + suffix).c_str()));
		return *this;
	}

	/**
	 * read a bamHeader to this file
	 * return the newly build header
	 */
	const BAMheader& readHeader() {
		return (header = BAMheader(sam_hdr_read(samFh)));
	}

	/**
	 * write the bamHeader to this file
	 */
	int writeHeader() const {
		return sam_hdr_write(samFh, header.bamHeader);
	}

	/**
	 * read next Bam record from a this SAMfile
	 */
	BAM read() {
		BAM aln;
		sam_read1(samFh, header.bamHeader, aln.bamAln);
		return aln;
	}

	/**
	 * write a BAM algnment to this SAMfile
	 */
	int write(const BAM& aln) {
		return sam_write1(samFh, header.bamHeader, aln.bamAln);
	}

	/* static methods */
	/**
	 * build and save an SAM/BAM index file for a given SAMfile basename, given suffix and given shift
	 */
	static int buildIndex(const string& basename, const string& suffix = "",
			int shift = DEFAULT_INDEX_SHIFT) {
		if(suffix.empty())
			return sam_index_build(basename.c_str(), shift);
		else
			return sam_index_build2(basename.c_str(), (basename + suffix).c_str(), shift);
	}

	/**
	 * build and save an SAM/BAM index file for a given SAMfile basename, given suffix, given shift and given threads
	 */
	static int buildIndex(const string& basename, int nThreads, const string& suffix = HTSindex::DEFAULT_INDEX_SUFFIX,
			int shift = DEFAULT_INDEX_SHIFT) {
		return sam_index_build3(basename.c_str(), (basename + suffix).c_str(), shift, nThreads);
	}

	/* static member fields */
	static const int DEFAULT_INDEX_SHIFT = 0;
	static const int DEFAULT_INDEX_BUILD_THREAD = 1;

	/* member fields */
private:
	samFile *samFh = nullptr;
	BAMheader header;
	HTSindex idx;
};

} /* namespace SAMtools */
} /* namespace EGriceLab */

#endif /* SRC_SAMFILE_H_ */
