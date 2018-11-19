/*
 * HTSindex.h
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */

#ifndef SRC_HTSINDEX_H_
#define SRC_HTSINDEX_H_

#include <string>
#include <cstdint>
#include <utility>
#include <htslib/hts.h>

namespace EGriceLab {
namespace SAMtools {
using std::string;

/**
 * A C++ wrapper class for the latest BAM/CAM index,
 * derived from with the __hts_idx_tsamFile struct described from htslib C library
 */
class HTSindex {
public:
	/* constructors */
	/** default constructor */
	HTSindex() = default;

	/** no copy constructor */
	HTSindex(const HTSindex&) = delete;
	/** no copy assignment */
	HTSindex& operator=(const HTSindex&) = delete;

	/** default move constructor */
	HTSindex(HTSindex&&) = default;

	/** move assignment */
	HTSindex& operator=(HTSindex&& other) {
		hts_idx_destroy(htsIndex);
		htsIndex = std::move(other.htsIndex);
		return *this;
	}

	/** destructor */
	virtual ~HTSindex() {
		hts_idx_destroy(htsIndex);
	}

	/** construct an HTSindex from given info */
	HTSindex(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls) {
		htsIndex = hts_idx_init(n, fmt, offset0, min_shift, n_lvls);
	}

	/** construct an HTSindex from a fresh copy of hts_idx_t */
	HTSindex(hts_idx_t* htsIndex) : htsIndex(htsIndex)
	{  	}

	/** construct a HTSindex from a file and format */
	explicit HTSindex(const string& basename, int fmt = DEFAULT_INDEX_FMT) {
		htsIndex = hts_idx_load(basename.c_str(), fmt);
	}

	/** construct a HTSindex from a file and suffix */
	HTSindex(const string& basename, const string& suffix) {
		htsIndex = hts_idx_load2(basename.c_str(), (basename + suffix).c_str());
	}

public:
	/* member methods */
	/** push a new content to the index */
    int push(int tid, int beg, int end, uint64_t offset, int is_mapped) {
    	return hts_idx_push(htsIndex, tid, beg, end, offset, is_mapped);
    }

    /** finish this index */
    void finish(uint64_t final_offset) {
    	hts_idx_finish(htsIndex, final_offset);
    }

    /** save this index to a given file, with optional index suffix */
    int save(const string& basename, const string& suffix = DEFAULT_INDEX_SUFFIX, int fmt = DEFAULT_INDEX_FMT) const {
    	if(suffix.empty())
    		return hts_idx_save(htsIndex, basename.c_str(), fmt);
    	else
    		return hts_idx_save_as(htsIndex, basename.c_str(), (basename + suffix).c_str(), fmt);
    }

    /* static methods */
    /** load a HTSindex from given file and format */
    static HTSindex load(const string& basename, int fmt = DEFAULT_INDEX_FMT) {
    	return HTSindex(basename, fmt);
    }

    /** load a HTSindex from given file and suffix */
    static HTSindex load(const string& basename, const string& suffix) {
    	return HTSindex(basename, suffix);
    }

    /* static fields */
	static const string DEFAULT_INDEX_SUFFIX; // default index format as .bai
	static const int DEFAULT_INDEX_FMT = 0; // 0 for .bai index

	/* member fields */
private:
	hts_idx_t *htsIndex = nullptr;

	friend class SAMfile;
};

} /* namespace SAMtools */
} /* namespace EGriceLab */

#endif /* SRC_HTSINDEX_H_ */
