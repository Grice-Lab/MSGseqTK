/*
 * BAM.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <cassert>
#include <boost/lexical_cast.hpp>
#include "BAM.h"

namespace EGriceLab {
namespace SAMtools {

BAM::BAM(const string& qname, uint16_t flag, int32_t tid, int32_t pos, uint8_t mapQ,
		const cigar_str& cigar, uint32_t l_seq, const seq_str& seq, const qual_str& qual,
		int32_t mtid, int32_t mpos, int32_t isize, uint64_t id) : BAM() {
	assert(seq.length() == (l_seq + 1) / 2);
	assert(qual.length() == l_seq);
	/* set core data */
	bamAln->core.l_qname = qname.length() + 1;
	bamAln->core.l_extranul = bamAln->core.l_qname % 4 != 0 ? 4 - bamAln->core.l_qname % 4 : 0;
	bamAln->core.l_qname += bamAln->core.l_extranul;
	bamAln->core.flag = flag;
	bamAln->core.pos = pos;
	bamAln->core.qual = mapQ;
	bamAln->core.n_cigar = cigar.length();
	bamAln->core.l_qseq = l_seq;
	bamAln->core.isize = isize;
	bamAln->core.mpos = mpos;
	bamAln->core.mtid = mtid;
#ifndef BAM_NO_ID
	bamAln->id = id;
#endif
	/* set variable length data and other fields */
	bamAln->l_data = bamAln->core.l_qname + bamAln->core.n_cigar * sizeof(uint32_t) + (bamAln->core.l_qseq + 1) / 2 + bamAln->core.l_qseq;
	bamAln->m_data = bamAln->l_data;
	bamAln->data = new uint8_t[bamAln->l_data]();
	uint8_t* ptr = bamAln->data;
	std::copy(qname.begin(), qname.end(), ptr); // copy qname
	ptr += bamAln->core.l_qname; // jump over qname and nulls
	std::copy_n(reinterpret_cast<const uint8_t*>(cigar.c_str()), bamAln->core.n_cigar * sizeof(uint32_t), ptr); // copy cigar
	ptr += bamAln->core.n_cigar * sizeof(uint32_t);
	std::copy(seq.begin(), seq.end(), ptr); // copy qseq
	ptr += seq.length();
	std::copy(qual.begin(), qual.end(), ptr); // copy qual
	ptr += qual.length();
	/* set bin */
	bamAln->core.bin = hts_reg2bin(getPos(), bamAln->core.pos + getAlignLen(), 14, 5); // 14, 5 are fixed values for BAM bins
}

void BAM::setQName(uint8_t l, const char* name) {
	uint8_t l_extranul = l % 4 != 0 ? 4 - l % 4 : 0; // extranulls needed for 32bit align
	int8_t l_diff = l + l_extranul - getQNameLen();
	uint32_t l_data = getDataLen() + l_diff;
	assert(l_data % 4 == 0);
	uint8_t* data = new uint8_t[l_data](); // value init
	uint8_t* ptr = data;
	std::copy(name, name + l, ptr); // copy new qname
	ptr += l + l_extranul; // jump over extra nulls
	std::copy(reinterpret_cast<const uint8_t*>(getCigar()), getData() + getDataLen(), ptr); // copy remaining old data
	delete[] bamAln->data; // delete old data
	bamAln->data = data; // use new data
	/* update length */
	setQNameLen(l + l_extranul);
	setExtranullLen(l_extranul);
	setDataLen(l_data);
}

void BAM::setCigar(uint32_t n, const uint32_t* cigar) {
	int32_t n_diff = n - getCigarNum();
	uint32_t l_data = getDataLen() + n_diff * sizeof(uint32_t);
	assert(l_data % 4 == 0);
	uint8_t* data = new uint8_t[l_data](); // value init
	uint8_t* ptr = data;
	std::copy(getData(), reinterpret_cast<const uint8_t*>(getCigar()), ptr); // copy old qname untill cigar
	ptr += getQNameLen();
	std::copy_n(reinterpret_cast<const uint8_t*>(cigar), n * sizeof(uint32_t), ptr); // copy new cigar
	ptr += n * sizeof(uint32_t);
	std::copy(getSeq(), getData() + getDataLen(), ptr); // copy remaining old data
	delete[] bamAln->data; // delete old data
	bamAln->data = data; // use new data
	/* update length */
	setCigarNum(n);
	setDataLen(l_data);
}

void BAM::setSeq(uint32_t l, const uint8_t* seq) {
	int32_t l_diff = l - getSeqLen();
	uint32_t l_data = getDataLen() + l_diff;
	assert(l_data % 4 == 0);
	uint8_t* data = new uint8_t[l_data];
	uint8_t* ptr = data;
	std::copy(getData(), getSeq(), ptr); // copy old data untill seq
	ptr += getQNameLen() + getCigarNum() * sizeof(uint32_t);
	std::copy_n(seq, l, ptr); // copy new seq
	ptr += l;
	std::copy(getQual(), getData() + getDataLen(), ptr); // copy remaining old data
	delete[] bamAln->data; // delete old data
	bamAln->data = data; // use new data
	/* update length */
	setSeqLen(l);
	setDataLen(l_data);
}

void BAM::setQual(uint32_t l, const uint8_t* qual) {
	if(l != getQualLen())
		return;
	int32_t l_diff = l - getQualLen();
	uint32_t l_data = getDataLen() + l_diff;
	assert(l_data % 4 == 0);
	uint8_t* data = new uint8_t[l_data];
	uint8_t* ptr = data;
	std::copy(getData(), getQual(), ptr); // copy old data untill qual
	ptr += getQNameLen() + getCigarNum() * sizeof(uint32_t) + (getSeqLen() + 1) / 2 /* ceil(seqLen / 2) */;
	std::copy_n(qual, l, ptr); // copy new qual
	ptr += getQualLen(); // jump over l_qseq
	std::copy(getAux(), getData() + getDataLen(), ptr); // copy remaining old data
	delete[] bamAln->data; // delete old data
	bamAln->data = data; // use new data
	/* update length */
	setDataLen(l_data);
}

BAM::seq_str BAM::nt16Encode(const string& seqRaw) {
	const uint32_t L = seqRaw.length();
	seq_str seq((L + 1) / 2, 0); // ceil(L / 2)
	for(uint32_t i = 0; i < L; i += 2)
		seq[i / 2] = nt16Encode(seqRaw[i]) << 4 | nt16Encode(seqRaw[i+1]); // seqRaw[i+1] always valid with the null terminal
	return seq;
}

string BAM::nt16Decode(const seq_str& seq) {
	const uint32_t L = seq.length();
	string seqRaw;
	seqRaw.reserve(L * 2);
	for(uint32_t i = 0; i < L * 2; ++i) {
		uint8_t b = getSeqBase(seq, i);
		if(i < L * 2 - 1 || b != 0) // only the last base of a seq is potentially null
			seqRaw.push_back(nt16Decode(b));
	}
	return seqRaw;
}

string BAM::decodeCigar(const cigar_str& cigar) {
	string str;
	for(cigar_str::value_type c : cigar)
		str += boost::lexical_cast<string>(bam_cigar_oplen(c)) + bam_cigar_opchr(c);
	return str;
}

} /* namespace SAMtools */
} /* namespace EGriceLab */
