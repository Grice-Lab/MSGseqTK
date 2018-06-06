/*
 * Loc.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: zhengqi
 */

#include "Loc.h"

namespace EGriceLab {
namespace MSGseqClean {

std::istream& operator>>(std::istream& in, Loc& loc) {
	in.ignore(1, '[');
	in >> loc.start;
	in.ignore(1, ',');
	in >> loc.end; /* automatical ignore additional space */
	in.ignore(1, ')');
	return in;
}

} /* namespace MSGseqClean */
} /* namespace EGriceLab */
