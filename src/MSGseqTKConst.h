/*
 * MSGseqCleanConst.h
 *
 *  Created on: May 16, 2018
 *      Author: zhengqi
 */

#ifndef SRC_MSGSEQTKCONST_H_
#define SRC_MSGSEQTKCONST_H_

#include <string>

namespace EGriceLab {
namespace MSGseqTK {

using std::string;

const string GZIP_FILE_SUFFIX = ".gz";
const string BZIP2_FILE_SUFFIX = ".bz2";
const string METAGENOME_FILE_SUFFIX = ".mtg";
const string FMINDEX_FILE_SUFFIX = ".fmidx";

const int MAX_NAME_LENGTH = 4096;
const double PHRED_LOG_BASE = 10;

}
}



#endif /* SRC_MSGSEQTKCONST_H_ */
