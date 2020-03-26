/*
 * BAMheader.h
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */

#ifndef SRC_BAMHEADER_H_
#define SRC_BAMHEADER_H_

#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cstdint>
#include <htslib/sam.h>
#include <htslib/hts.h>

namespace EGriceLab {
namespace SAMtools {
using std::string;
using std::vector;

/**
 * A C++ wrapper class for the alignment header,
 * interfacing with the bam_hdr_t struct described from htslib C library
 */
class BAMheader {
public:
	typedef std::map<string, uint32_t> targetMap; /* target/chromosome length map */
	typedef std::map<string, string> textMap;     /* other plain-text map */
	/* constructors */
	/** default constructor */
	BAMheader() : bamHeader(bam_hdr_init())
	{  	}

	/** copy constructor */
	BAMheader(const BAMheader& other) : bamHeader(bam_hdr_dup(other.bamHeader))
	{  	}

	/** copy assignment */
	BAMheader& operator=(const BAMheader& other) {
		bam_hdr_destroy(bamHeader);
		bamHeader = bam_hdr_dup(other.bamHeader);
		return *this;
	}

	/** default move constructor */
	BAMheader(BAMheader&&) = default;

	/** move assignment */
	BAMheader& operator=(BAMheader&& other) {
		bam_hdr_destroy(bamHeader);
		bamHeader = std::move(other.bamHeader);
		return *this;
	}

	/** destructor */
	virtual ~BAMheader() {
		bam_hdr_destroy(bamHeader);
	}

	/** construct a BAMheader from a given bam_hdr_t */
	BAMheader(bam_hdr_t* bamHeader) : bamHeader(bamHeader)
	{  	}

	/** construct a BAMheader from a target map */
	BAMheader(const targetMap& targetDict) : BAMheader() {
		setTarget(targetDict);
	}

	/** construct a BAMheader from a given order of targets and a target map */
	BAMheader(const vector<string>& targets, const targetMap& targetDict) : BAMheader() {
		setTarget(targets, targetDict);
	}

	/* member methods */
	/** get index of a given chrom */
	int getIndex(const string& name) {
		return bam_name2id(bamHeader, name.c_str());
	}

	/** set/add the @HD tag value from this BAMheader */
	int setHDTag(const string& tag, const string& val) const {
		return sam_hdr_change_HD(bamHeader, tag.c_str(), val.c_str());
	}

	/** add a new target with given name and length */
	BAMheader& addTarget(const string& name, uint32_t len);

	/** set targets using a target dict */
	BAMheader& setTarget(const targetMap& targetDict);

	/** set targets using a given order of targets and a target dict */
	BAMheader& setTarget(const vector<string>& targets, const targetMap& targetDict);

	/** add a new text tag into this header */
	BAMheader& addTag(const string& tag, const string& val);

	/** set text tags using a text dict */
	BAMheader& setTag(const textMap& textDict);

	/** set text tags using a given order of tags and a text dict */
	BAMheader& setTag(const vector<string>& tags, const textMap& textDict);

	/* member fields */
private:
	bam_hdr_t *bamHeader = nullptr;

	friend class SAMfile;

	/* static fields */
public:
	static const char TEXT_TAG_SEP = '\t';
};

} /* namespace SAMtools */
} /* namespace EGriceLab */

#endif /* SRC_BAMHEADER_H_ */
