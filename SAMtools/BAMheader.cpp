/*
 * BamHeader.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */

#include <string>
#include <algorithm>
#include <utility>
#include "htslib/khash.h"
#include "BAMheader.h"

namespace EGriceLab {
namespace SAMtools {
using std::string;

BAMheader::BAMheader(const targetMap& targetDict) : BAMheader() {
	/* init targets */
	bamHeader->n_targets = targetDict.size();
	bamHeader->target_len = new uint32_t[bamHeader->n_targets];
	bamHeader->target_name = new char*[bamHeader->n_targets];
	/* copy each targets */
	uint32_t i = 0;
	for(const targetMap::value_type& target : targetDict) {
		const string& tname = target.first;
		const uint32_t tlen = target.second;
		bamHeader->target_len[i] = tlen;
		bamHeader->target_name[i] = new char[tname.length()];
		std::copy(tname.begin(), tname.end(), bamHeader->target_name[i]);
		i++;
	}
}

BAMheader::BAMheader(const targetMap& targetDict, const textMap& textDict) : BAMheader(targetDict) {
	/* calculate total text length */
	bamHeader->l_text = 0;
	for(const textMap::value_type& text : textDict)
		bamHeader->l_text += text.first.length() + text.second.length() + 2; /* with a sep and a new line */
	bamHeader->text = new char[bamHeader->l_text];
	char* ptr = bamHeader->text;
	for(const textMap::value_type& text : textDict) {
		std::copy(text.first.begin(), text.first.end(), ptr);
		ptr += text.first.length();
		*ptr++ = TEXT_TAG_SEP;
		std::copy(text.second.begin(), text.second.end(), ptr);
		ptr += text.second.length();
		*ptr++ = '\n';
	}
}

BAMheader& BAMheader::addTag(const string& tag, const string& val) {
	uint32_t l_text = bamHeader->l_text + tag.length() + val.length() + 2; // new l_text
	char* text = new char[l_text]; // next text
	char* ptr = text;
	std::copy(bamHeader->text, bamHeader->text + bamHeader->l_text, ptr); // copy old text
	ptr += bamHeader->l_text;

	std::copy(tag.begin(), tag.end(), ptr);
	ptr += tag.length();
	*ptr++ = TEXT_TAG_SEP;
	std::copy(val.begin(), val.end(), ptr);
	ptr += val.length();
	*ptr++ = '\n';

	/* repace data */
	bamHeader->l_text = l_text;
	delete[] bamHeader->text;
	bamHeader->text = text;
	return *this;
}

} /* namespace SAMtools */
} /* namespace EGriceLab */
