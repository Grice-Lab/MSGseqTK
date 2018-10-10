/*
 * BitSeqRRR.h
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#ifndef BITSEQRRR_H_
#define BITSEQRRR_H_

#include "BitSeq.h"
#include "BitStr.h"

namespace EGriceLab {
namespace libSDS {

/**
 * A Raman, Raman and Rao's implementation of BitSeq, based on the algorithm described in
 *  [1] R. Raman, V. Raman and S. Rao. Succinct indexable dictionaries with applications
 *     to encoding $k$-ary trees and multisets. SODA02.
 *  [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 */
class BitSeqRRR: public BitSeq {
public:
	/* nested types and enums */
	/**
	 * pre-computed table offsets used in BitSequenceRRR algorithm
	 */
	class TableOffset {

	};

	BitSeqRRR() {
		// TODO Auto-generated constructor stub

	}
	virtual ~BitSeqRRR() {
		// TODO Auto-generated destructor stub
	}

	/* member fields */
private:


};

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSEQRRR_H_ */
