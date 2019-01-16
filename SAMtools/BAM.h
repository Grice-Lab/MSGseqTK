/*
 * BAM.h
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */

#ifndef SRC_BAM_H_
#define SRC_BAM_H_

#include <string>
#include <utility>
#include <boost/algorithm/string/regex.hpp>
#include <htslib/sam.h>

namespace EGriceLab {
namespace SAMtools {
using std::string;
using std::basic_string;

/**
 * A C++ wrapper class for one bam alignment,
 * derived from the bam1_t struct described from htslib C library
 */
class BAM {
public:
	typedef string qname_str; /* qname is just a string */
	typedef basic_string<uint32_t> cigar_str; /* wrapper for encoded cigar str */
	typedef basic_string<uint8_t> seq_str; /* wrapper for 4-bits encode seq str */
	typedef basic_string<uint8_t> qual_str; /* wrapper for mapping quality str */

	/* constructors */
	/** default constructor */
	BAM() : bamAln(bam_init1())
	{  	}

	/** copy constructor */
	BAM(const BAM& other) {
		bam_copy1(bamAln, other.bamAln);
	}

	/** copy assignment */
	BAM& operator=(const BAM& other) {
		bam_destroy1(bamAln);
		bamAln = bam_init1();
		bam_copy1(bamAln, other.bamAln);
		return *this;
	}

	/** move constructor */
	BAM(BAM&&) = default;

	/** move assignment */
	BAM& operator=(BAM&& other) {
		bam_destroy1(bamAln);
		bamAln = std::move(other.bamAln);
		return *this;
	}

	/** construct a BAM from all core data */
	BAM(const string& qname, uint16_t flag, int32_t tid, int32_t pos, uint8_t mapQ,
			const cigar_str& cigar, uint32_t l_seq, const seq_str& seq, const qual_str& qual,
			int32_t mtid = -1, int32_t mpos = -1, int32_t isize = 0, uint64_t id = 0);

	/** deligating construt a BAM from all core data, using raw seq */
	BAM(const string& qname, uint16_t flag, int32_t tid, int32_t pos, uint8_t mapQ,
			const cigar_str& cigar, const string& seq, const qual_str& qual,
			int32_t mtid = -1, int32_t mpos = -1, int32_t isize = 0, uint64_t id = 0)
	: BAM(qname, flag, tid, pos, mapQ, cigar, seq.length(), nt16Encode(seq), qual, mtid, mpos, isize, id)
	{  }

	/* member methods */
	/** getters and setters for bam1_core_t member fields */
	/** get target/reference ID */
	int32_t getRId() const {
		return bamAln->core.tid;
	}

	/** set target/reference ID */
	void setRId(int32_t tid) {
		bamAln->core.tid = tid;
	}

	/** get pos 0-based */
	int32_t getPos() const {
		return bamAln->core.pos;
	}

	/** set pos */
	void setPos(int32_t pos) {
		bamAln->core.pos = pos;
	}

	/** get bin */
	uint16_t getBin() const {
		return bamAln->core.bin;
	}

	/** set bin */
	void setBin(uint16_t bin) {
		bamAln->core.bin = bin;
	}

	/** get qname len, including tailing nulls */
	uint8_t getQNameLen() const {
		return bamAln->core.l_qname;
	}

	/** set qname len */
	void setQNameLen(uint8_t len) {
		bamAln->core.l_qname = len;
	}

	/** get number of extra nulls */
	uint8_t getExtranulLen() const {
		return bamAln->core.l_extranul;
	}

	/** set number of extra nulls */
	void setExtranullLen(uint8_t len) {
		bamAln->core.l_extranul = len;
	}

	/** get number of cigar elements */
	uint32_t getCigarNum() const {
		return bamAln->core.n_cigar;
	}

	/** set number of cigar elements */
	void setCigarNum(uint32_t n_cigar) {
		bamAln->core.n_cigar = n_cigar;
	}

	/** get qseq length */
	int32_t getSeqLen() const {
		return bamAln->core.l_qseq;
	}

	/** set qseq length */
	void setSeqLen(int32_t l_qseq) {
		bamAln->core.l_qseq = l_qseq;
	}

	/** get total data len */
	uint32_t getDataLen() const {
		return bamAln->l_data;
	}

	/** set total data len */
	void setDataLen(uint32_t len) {
		bamAln->l_data = len;
	}

	/** get qual len, alias of getSeqLen() */
	int32_t getQualLen() const {
		return getSeqLen();
	}

	/** get flag */
	uint16_t getFlag() const {
		return bamAln->core.flag;
	}

	/** set flag */
	void setFlag(uint16_t flag) {
		bamAln->core.flag = flag;
	}

	bool getPairedFlag() const {
		return bamAln->core.flag & BAM_FPAIRED != 0;
	}

	void setPairedFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FPAIRED;
		else
			bamAln->core.flag &= ~BAM_FPAIRED;
	}

	bool getProperPairFlag() const {
		return bamAln->core.flag & BAM_FPROPER_PAIR != 0;
	}

	void setProperPairFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FPAIRED;
		else
			bamAln->core.flag &= ~BAM_FPAIRED;
	}

	bool getUnmapFlag() const {
		return bamAln->core.flag & BAM_FUNMAP != 0;
	}

	void setUnmapFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FUNMAP;
		else
			bamAln->core.flag &= ~BAM_FUNMAP;
	}

	bool getMateUnmapFlag() const {
		return bamAln->core.flag & BAM_FMUNMAP != 0;
	}

	void setMateUnmapFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FMUNMAP;
		else
			bamAln->core.flag &= ~BAM_FMUNMAP;
	}

	bool getRevcomFlag() const {
		return bamAln->core.flag & BAM_FREVERSE != 0;
	}

	void setRevcomFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FREVERSE;
		else
			bamAln->core.flag &= ~BAM_FREVERSE;
	}

	bool getMateReverseFlag() const {
		return bamAln->core.flag & BAM_FMREVERSE != 0;
	}

	void setMateRevcomFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FMREVERSE;
		else
			bamAln->core.flag &= ~BAM_FMREVERSE;
	}

	bool getIsFwdFlag() const {
		return bamAln->core.flag & BAM_FREAD1 != 0;
	}

	void setIsFwdFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FREAD1;
		else
			bamAln->core.flag &= ~BAM_FREAD1;
	}

	bool getIsRevFlag() const {
		return bamAln->core.flag & BAM_FREAD2 != 0;
	}

	void setIsRevFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FREAD2;
		else
			bamAln->core.flag &= ~BAM_FREAD2;
	}

	bool getSecondaryFlag() const {
		return bamAln->core.flag & BAM_FSECONDARY != 0;
	}

	void setSecondaryFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FSECONDARY;
		else
			bamAln->core.flag &= ~BAM_FSECONDARY;
	}

	bool getQCFailFlag() const {
		return bamAln->core.flag & BAM_FQCFAIL != 0;
	}

	void setQCFailFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FQCFAIL;
		else
			bamAln->core.flag &= ~BAM_FQCFAIL;
	}

	bool getDupFlagFlag() const {
		return bamAln->core.flag & BAM_FDUP != 0;
	}

	void setDupFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FDUP;
		else
			bamAln->core.flag &= ~BAM_FDUP;
	}

	bool getSupplementaryFlag() const {
		return bamAln->core.flag & BAM_FSUPPLEMENTARY != 0;
	}

	void setSupplementaryFlag(bool flag = true) {
		if(flag)
			bamAln->core.flag |= BAM_FSUPPLEMENTARY;
		else
			bamAln->core.flag &= ~BAM_FSUPPLEMENTARY;
	}

	/** methods for other fields */
#ifndef BAM_NO_ID
	/** get align ID */
	uint64_t getId() const {
		return bamAln->id;
	}

	/** set align id */
	void setId(uint64_t id) {
		bamAln->id = id;
	}
#endif

	/** get the entire data */
	const uint8_t* getData() const {
		return bamAln->data;
	}

	/** Get the qname */
	const char* getQName() const {
		return bam_get_qname(bamAln);
	}

	/** get the qname_str */
	qname_str getQNameStr() const {
		return qname_str(bam_get_qname(bamAln), bamAln->core.l_qname - bamAln->core.l_extranul); // l_qname includes tailing nulls
	}

	/**
	 * set qname, need to ducplicate the old data and makes it 32bit aligned
	 * @param len  length of new name, !! includes the null terminal !!
	 * @param name  new name_str
	 */
	void setQName(uint8_t l, const char* name);

	/** set qnameStr */
	void setQName(const qname_str& name) {
		setQName(name.length() + 1, name.c_str());
	}

	/**
	 Get whether the query is on the reverse strand
	 @param  b  pointer to an alignment
	 @return    true if query is on the reverse strand
	 */
	bool isRev() const {
		return bam_is_rev(bamAln);
	}

	/**
 	 Get whether the query's mate is on the reverse strand
 	 b  pointer to an alignment
 	 @return    true if query's mate on the reverse strand
	 */
	bool isMateRev() const {
		return bam_is_mrev(bamAln);
	}

	/**
	 Get the CIGAR array
	 @return    pointer to the CIGAR array
	 In the CIGAR array, each element is a 32-bit integer. The
	 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
	 length of a CIGAR.
	 */
	const uint32_t* getCigar() const {
		return bam_get_cigar(bamAln);
	}

	/** get the cigar_str */
	cigar_str getCigarStr() const {
		return cigar_str(bam_get_cigar(bamAln), bamAln->core.n_cigar);
	}

	/** get the decoded cigar string */
	string getCigarStrDecoded() const {
		return decodeCigar(getCigarStr());
	}

	/**
	 * get a cigar op at given pos
	 * @return  a uint32_t encoded cigarOp
	 */
	uint32_t getCigarOp(uint32_t i) const {
		return bam_cigar_op(bam_get_cigar(bamAln)[i]);
	}

	/**
	 * get a cigar op len at given pos
	 * @return  a the cigarOp len
	 */
	uint32_t getCigarOpLen(uint32_t i) const {
		return bam_cigar_oplen(bam_get_cigar(bamAln)[i]);
	}

	/**
	 * get a cigar op at given pos
	 * @return  a decoded cigarOp character, one of "MIDNSHP=XB", or '?' if failed
	 */
	char getCigarOpChar(uint32_t i) const {
		return bam_cigar_opchr(bam_get_cigar(bamAln)[i]);
	}

	/**
	 * set the cigar, need to ducplicate the entire data
	 * @param n  elements of new cigar, encoded as uint32_t arrays
	 * @param cigar  new cigar array
	 */
	void setCigar(uint32_t n, const uint32_t* cigar);

	/**
	 * set the cigar, need to ducplicate the entire data
	 * @param cigar  new cigar basic_string
	 */
	void setCigar(const cigar_str& cigar) {
		setCigar(cigar.length(), cigar.c_str());
	}

	/**
	 * set the cigar string, need to ducplicate the entire data
	 * @param cigar  an encoded cigar string
	 */
	void setCigarStr(const string& cigarStr) {
		setCigar(encodeCigar(cigarStr));
	}

	/**
	 * Get query sequence
	 * @return    pointer to sequence
	 * Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
	 * 8 for T and 15 for N. Two bases are packed in one byte with the base
	 * at the higher 4 bits having smaller coordinate on the read. It is
	 * recommended to use bam_seqi() macro to get the base.
	 */
	const uint8_t* getSeq() const {
		return bam_get_seq(bamAln);
	}

	/**
	 * get the encoded seq_str
	 */
	const seq_str getSeqStr() const {
		return seq_str(bam_get_seq(bamAln), bamAln->core.l_qseq);
	}


	/**
	 * get the decoded seq string
	 */
	const string getSeqStrDecoded() const {
		return nt16Decode(getSeqStr());
	}

	/**
	 * set query sequence
	 * @param n  length of the encoded seq
	 * @param seq  the new encoded seq
	 */
	void setSeq(uint32_t l, const uint8_t* seq);

	/**
	 * set query sequence
	 * @param seq  the new encoded seq_str
	 */
	void setSeq(const seq_str& seq) {
		setSeq(seq.length(), seq.c_str());
	}

	/**
	 * set query sequence
	 * @param seqStr  un-encoded seq
	 */
	void setSeq(const string& seq) {
		setSeq(nt16Encode(seq));
	}

	/**
	 * Get query quality
	 * @return    pointer to quality string
	 */
	const uint8_t* getQual() const {
		return bam_get_qual(bamAln);
	}

	/** get quality str */
	qual_str getQualStr() const {
		return qual_str(bam_get_qual(bamAln), bamAln->core.l_qseq); // l_qseq is l_qual
	}

	/** set qual */
	void setQual(uint32_t l, const uint8_t* qual);

	/** set qual with qual_str */
	void setQual(const qual_str& qual) {
		setQual(qual.length(), qual.c_str());
	}

	/**
	 * Get auxiliary data
	 * @return    pointer to the concatenated auxiliary data,
	 * which has the structure: qname-cigar-seq-qual-aux
	 */
	const uint8_t* getAux() const {
		return bam_get_aux(bamAln);
	}

	/** get length of auxiliary data */
	uint32_t getAuxLen() const {
		return bam_get_l_aux(bamAln);
	}

	/** get a aux of int type */
	int64_t getAuxInt(const string& tag) const {
		return bam_aux2i(bam_aux_get(bamAln, tag.c_str()));
	}

	/** get a aux of float type */
	double getAuxFloat(const string& tag) const {
		return bam_aux2f(bam_aux_get(bamAln, tag.c_str()));
	}

	/** get a aux of char type */
	char getAuxChar(const string& tag) const {
		return bam_aux2A(bam_aux_get(bamAln, tag.c_str()));
	}

	/** get a aux of string type */
	string getAuxStr(const string& tag) const {
		return bam_aux2Z(bam_aux_get(bamAln, tag.c_str()));
	}

	/** get a aux of array of int type at position i */
	int64_t getAuxInt(const string& tag, uint32_t i) const {
		return bam_auxB2i(bam_aux_get(bamAln, tag.c_str()), i);
	}

	/** get a aux of array of float type at position i */
	double getAuxFloat(const string& tag, uint32_t i) const {
		return bam_auxB2f(bam_aux_get(bamAln, tag.c_str()), i);
	}

	/** add or update a new int tag for this BAM record */
	int setAux(const string& tag, int32_t val) {
		return bam_aux_update_int(bamAln, tag.c_str(), val);
	}

	/** add or update a new int tag for this BAM record */
	int setAux(const string& tag, uint32_t val) {
		return bam_aux_update_int(bamAln, tag.c_str(), val);
	}

	/** add or update a new int tag for this BAM record */
	int setAux(const string& tag, int64_t val) {
		return bam_aux_update_int(bamAln, tag.c_str(), val);
	}

	/** add or update a new int tag for this BAM record */
	int setAux(const string& tag, uint64_t val) {
		return bam_aux_update_int(bamAln, tag.c_str(), val);
	}

	/** add or update a new float tag for this BAM record */
	int setAux(const string& tag, float val) {
		return bam_aux_update_float(bamAln, tag.c_str(), val);
	}

	/** add or update a new float tag for this BAM record */
	int setAux(const string& tag, double val) {
		return bam_aux_update_float(bamAln, tag.c_str(), static_cast<float>(val));
	}

	/** add or update a new string tag for this BAM record */
	int setAux(const string& tag, const string& val) {
		return bam_aux_update_str(bamAln, tag.c_str(), val.length() + 1, val.c_str()); // null included
	}

	/** add or update a new array type aux record for this BAM */
	template<typename T>
	int setAux(const string& tag, uint8_t type, uint32_t nItems, T* data) {
		return bam_aux_update_array(bamAln, tag.c_str(), type, nItems, (void *) data);
	}

	/** remove a aux tag from this BAM */
	int removeAux(const string& tag) {
		return bam_aux_del(bamAln, bam_aux_get(bamAln, tag.c_str()));
	}

	/** get query length from cigar */
	int getCigarQLen() const {
		return bam_cigar2qlen(bamAln->core.n_cigar, getCigar());
	}

	/** get reference length from cigar */
	int getCigarRLen() const {
		return bam_cigar2rlen(bamAln->core.n_cigar, getCigar());
	}

	/** get alignment length, alias to getCigarRLen() */
	int getAlignLen() const {
		return getCigarRLen();
	}

	/** get alignment end, 1-based */
	int32_t getAlignEnd() const {
		return getPos() + getAlignLen();
	}

	/* member fields */
private:
	bam1_t *bamAln = nullptr;

public:
	friend class SAMfile;

	/* static methods */
	/**
	 * converting a nucleotide character to 4-bit encoding.
	 * The input character may be either an IUPAC ambiguity code, '=' for 0, or
	 * '0'/'1'/'2'/'3' for a result of 1/2/4/8.  The result is encoded as 1/2/4/8
	 * for A/C/G/T or combinations of these bits for ambiguous bases.
	 */
	static uint8_t nt16Encode(char c) {
		return seq_nt16_table[c];
	}

	/**
	 * converting a 4-bit encoded nucleotide to an IUPAC
	 * ambiguity code letter (or '=' when given 0).
	 */
	static char nt16Decode(uint8_t b) {
		return seq_nt16_str[b];
	}

	/** get the i-th base of an encoded seq */
	static uint8_t getSeqBase(uint8_t* seq, uint32_t i) {
		return bam_seqi(seq, i);
	}

	/** get the i-th base of an encoded seq_str */
	static uint8_t getSeqBase(const seq_str& seq, uint32_t i) {
		return bam_seqi(seq, i);
	}

	/**
	 * encode raw seq to a string
	 */
	static seq_str nt16Encode(const string& seqRaw);

	/**
	 * decode encoded seq_str to raw string with given length
	 */
	static string nt16Decode(const seq_str& seq);

	/** encode cigar operator from character to int
	 * @return -1  if not valid
	 */
	static int encodeCigar(char op) {
		return string(BAM_CIGAR_STR).find(op);
	}

	/*
	 * encode readable cigar string to cigar_str
	 * @param  str  readable cigar string
	 * @return  encoded cigar_str
	 */
	static cigar_str encodeCigar(const string& str);

	/*
	 * encode readable cigar string to cigar_str
	 * @param  cigar  encoded cigar_str
	 * @return  decoded readable cigar string
	 */
	static string decodeCigar(const cigar_str& cigar);

	/** get a cigar type
	 * @param op  cigar operator
	 * @return a bitmask with bit 1: consum query; bit 2: consume reference
	 */
	static uint32_t getCigarType(uint32_t op) {
		return bam_cigar_type(op);
	}

	/** get query length from a given cigar */
	static int cigar2QLen(uint32_t n, const uint32_t* cigar) {
		return bam_cigar2qlen(n, cigar);
	}

	/** get query length from a given cigar_str */
	static int cigar2QLen(const cigar_str& cigar) {
		return cigar2QLen(cigar.length(), cigar.c_str());
	}

	/** get reference length from cigar */
	static int cigar2RLen(uint32_t n, const uint32_t* cigar) {
		return bam_cigar2rlen(n, cigar);
	}

	/** get reference length from cigar_str */
	static int cigar2RLen(const cigar_str& cigar) {
		return cigar2RLen(cigar.length(), cigar.c_str());
	}

	/* static fields */
	static const boost::regex CIGAR_PATTERN;
	static const double DEFAULT_AUX_DATA_FACTOR; // data factor between m_data and l_data, to avoid too much internal re-allocation by htslib
};

} /* namespace SAMtools */
} /* namespace EGriceLab */

#endif /* SRC_BAM_H_ */
