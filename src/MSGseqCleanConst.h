/*
 * MSGseqCleanConst.h
 *
 *  Created on: May 16, 2018
 *      Author: zhengqi
 */

#ifndef SRC_MSGSEQCLEANCONST_H_
#define SRC_MSGSEQCLEANCONST_H_

#include <string>

namespace EGriceLab {
namespace MSGseqClean {

using std::string;

const string GZIP_FILE_SUFFIX = ".gz";
const string BZIP2_FILE_SUFFIX = ".bz2";
const string METAGENOME_FILE_SUFFIX = ".mtg";
const string RFMINDEX_FILE_SUFFIX = ".rfm";

const int MAX_NAME_LENGTH = 4096;

}
}



#endif /* SRC_MSGSEQCLEANCONST_H_ */
