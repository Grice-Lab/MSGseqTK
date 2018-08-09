/*
 * GFF3.h
 *
 *  Created on: Aug 9, 2018
 *      Author: zhengqi
 */

#ifndef UCSC_GFF3_H_
#define UCSC_GFF3_H_

#include "GFF.h"

namespace EGriceLab {
namespace UCSC {

class GFF3: public GFF {
public:
	/* constructors */
	/** default constructor */
	GFF3() {  }

	/** virtual destructor */
	virtual ~GFF3() {  }

	/* implementation of base class abstract methods */
	/** read attrs from formated text */
	virtual void readAttributes(const string& attrStr);

	/** write attrs to formated text */
	virtual string writeAttributes() const;
};

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* UCSC_GFF3_H_ */
