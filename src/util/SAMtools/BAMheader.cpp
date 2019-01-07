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

BAMheader& BAMheader::addTarget(const string& name, uint32_t len) {
	/* init new target values */
	uint32_t n_targets = bamHeader->n_targets + 1;
	uint32_t* target_len = new uint32_t[n_targets];
	char **target_name = new char*[n_targets];

	/* copy old values */
	std::copy_n(bamHeader->target_len, n_targets, target_len);
	std::copy_n(bamHeader->target_name, n_targets, target_name);

	/* add new value */
	target_len[n_targets - 1] = len;
	target_name[n_targets - 1] = new char[name.length() + 1]; // include null-terminal
	std::copy_n(name.c_str(), name.length() + 1, target_name[n_targets - 1]);

	/* replace data */
	bamHeader->n_targets = n_targets;
	delete[] bamHeader->target_len;
	bamHeader->target_len = target_len;
	delete[] bamHeader->target_name;
	bamHeader->target_name = target_name;

	return *this;
}

BAMheader& BAMheader::setTarget(const targetMap& targetDict) {
	/* delete old data */
	delete[] bamHeader->target_len;
	delete[] bamHeader->target_name;

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
		bamHeader->target_name[i] = new char[tname.length() + 1]; /* including null in target name */
		std::copy_n(tname.c_str(), tname.length() + 1, bamHeader->target_name[i]);
		i++;
	}

	return *this;
}

BAMheader& BAMheader::setTarget(const vector<string>& targets, const targetMap& targetDict) {
	/* delete old data */
	delete[] bamHeader->target_len;
	delete[] bamHeader->target_name;

	/* init targets */
	bamHeader->n_targets = targets.size();
	bamHeader->target_len = new uint32_t[bamHeader->n_targets];
	bamHeader->target_name = new char*[bamHeader->n_targets];
	/* copy each targets */
	for(uint32_t i = 0; i < targets.size(); ++i) {
		const string& tname = targets[i];
		const uint32_t tlen = targetDict.at(tname);
		bamHeader->target_len[i] = tlen;
		bamHeader->target_name[i] = new char[tname.length() + 1]; /* including null in target name */
		std::copy_n(tname.c_str(), tname.length() + 1, bamHeader->target_name[i]);
	}

	return *this;
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

BAMheader& BAMheader::setTag(const textMap& textDict) {
	/* delete old data */
	delete[] bamHeader->text;

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

	return *this;
}

BAMheader& BAMheader::setTag(const vector<string>& tags, const textMap& textDict) {
	/* delete old data */
	delete[] bamHeader->text;

	/* calculate total text length */
	bamHeader->l_text = 0;
	for(const string& tag : tags)
		bamHeader->l_text += tag.length() + textDict.at(tag).length() + 2; /* with a sep and a new line */
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

	return *this;
}

} /* namespace SAMtools */
} /* namespace EGriceLab */
