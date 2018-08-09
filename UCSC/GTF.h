/*
 * GTF.h
 *
 *  Created on: Aug 9, 2018
 *      Author: zhengqi
 */

#ifndef UCSC_GTF_H_
#define UCSC_GTF_H_

#include "GFF.h"

namespace EGriceLab {
namespace UCSC {

class GTF: public GFF {
public:
	/* constructors */
	/** default constructor */
	GTF() {  }

	/** virtual destructor */
	virtual ~GTF() {  }

	/* implementation of base class abstract methods */
	/** read attrs from formated text */
	virtual void readAttributes(const string& attrStr);

	/** write attrs to formated text */
	virtual string writeAttributes() const;
};

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* UCSC_GTF_H_ */
