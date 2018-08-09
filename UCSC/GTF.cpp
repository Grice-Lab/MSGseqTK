/*
 * GTF.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: zhengqi
 */

#include <boost/algorithm/string.hpp>
#include "GTF.h"

namespace EGriceLab {
namespace UCSC {
using namespace std;

void GTF::readAttributes(const string& attrStr) {
	vector<string> attrs;
	boost::split(attrs, attrStr, boost::is_any_of(" \";"), boost::token_compress_on);
	for(vector<string>::size_type i = 0; i < attrs.size(); i += 2)
		setAttr(attrs[i], attrs[i+1]);
}

string GTF::writeAttributes() const {
	string attrStr;
	const attr_map& attrValues = getAttrValues();
	for(attr_map::const_iterator pair = attrValues.begin(); pair != attrValues.end(); ++pair) {
		if(pair != attrValues.begin()) /* non-first */
			attrStr += "; ";
		attrStr += pair->first + " \"" + pair->second + "\"";
	}
	return attrStr;
}

} /* namespace UCSC */
} /* namespace EGriceLab */


