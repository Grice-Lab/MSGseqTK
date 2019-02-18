/*
 * MSGseqTKConst.h
 *
 *  Created on: May 16, 2018
 *      Author: zhengqi
 */

#ifndef SRC_MSGSEQTKCONST_H_
#define SRC_MSGSEQTKCONST_H_

#include <string>
#include <limits>
#include <cfloat>
#include <cassert>

using std::string;

namespace EGriceLab {
namespace MSGseqTK {

const double inf = std::numeric_limits<double>::infinity();
const double infV = -inf;
const std::streamsize MAX_STREAM_SIZE = std::numeric_limits<std::streamsize>::max();

const string GZIP_FILE_SUFFIX = ".gz";
const string BZIP2_FILE_SUFFIX = ".bz2";
const string METAGENOME_FILE_SUFFIX = ".mtg";
const string FMINDEX_FILE_SUFFIX = ".fmidx";
const string FMDINDEX_FILE_SUFFIX = ".fmdi";
const string SAM_SUFFIX = ".sam";
const string BAM_SUFFIX = ".bam";

const int MAX_NAME_LENGTH = 4096;
const double PHRED_LOG_BASE = 10;

const double MIN_LOGLIK_EXP = DBL_MIN_EXP; // min loglik exp allowed to avoid numeric under-flow

}
}



#endif /* SRC_MSGSEQTKCONST_H_ */
