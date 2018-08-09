/*
 * GFF.h
 *
 *  Created on: Aug 7, 2018
 *      Author: zhengqi
 */

#ifndef GFF_H_
#define GFF_H_
#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace EGriceLab {
namespace UCSC {

using std::string;
using std::vector;
using std::map;
using std::istream;
using std::ostream;
using std::cerr;
using std::endl;

/**
 * an abstract class represent a UCSC GFF (GTF or GFF3) record
 */
class GFF {
public:
	typedef vector<string> attr_list;
	typedef map<string, string> attr_map;

	enum VERSION { UNK = 1, GTF, GFF3 };

	/** virtual destructor */
	virtual ~GFF() {  }

	const attr_map& getAttrValues() const {
		return attrValues;
	}

	long getEnd() const {
		return end;
	}

	void setEnd(long end) {
		this->end = end;
	}

	int getFrame() const {
		return frame;
	}

	void setFrame(int frame) {
		this->frame = frame;
	}

	double getScore() const {
		return score;
	}

	void setScore(double score) {
		this->score = score;
	}

	const string& getSeqname() const {
		return seqname;
	}

	void setSeqname(const string& seqname) {
		this->seqname = seqname;
	}

	const string& getSource() const {
		return source;
	}

	void setSource(const string& source) {
		this->source = source;
	}

	long getStart() const {
		return start;
	}

	void setStart(long start) {
		this->start = start;
	}

	char getStrand() const {
		return strand;
	}

	void setStrand(char strand) {
		this->strand = strand;
	}

	const string& getType() const {
		return type;
	}

	void setType(const string& type) {
		this->type = type;
	}

	/* member methods */
	/** get number of attrs */
	size_t numAttrs() const {
		return attrValues.size();
	}

	/** get the attr value given its name */
	string getAttr(const string& name) const {
		return attrValues.at(name);
	}

	/** set the attr value with given name, create if not already exists */
	void setAttr(const string& name, const string& val) {
		attrValues[name] = val;
	}

	/** test whether this attr name exists */
	bool hasAttr(const string& name) const {
		return attrValues.count(name) > 0;
	}

	/** load object from binary input */
	istream& load(istream& in);

	/** save object to binary output */
	ostream& save(ostream& out) const;

	/* abstract methods */
	/** read attrs from formated text */
	virtual void readAttributes(const string& attrStr) = 0;

	/** write attrs to formated text */
	virtual string writeAttributes() const = 0;

	/* non-member operators */
	/* unformated read */
	friend istream& operator>>(istream& in, GFF& record);

	/* unformated write */
	friend ostream& operator<<(ostream& out, const GFF& record);

	/* relationship operators */
	/** compare whether two GFF objects is the same
	 * @return true  if and only if all fields except the attrs are equal
	 */
	friend bool operator==(const GFF& lhs, const GFF& rhs);

private:
	/* member fields */
	string seqname;
	string source;
	string type;
	long start;
	long end;
	double score;
	char strand;
	int frame;
	attr_map attrValues; /* attribute name->value map */

public:
	/* class constants */
	static const char SEP = '\t';
	static const string GFF_SUFFIX;
	static const string GTF_SUFFIX;
	static const string GFF3_SUFFIX;
	static const double INVALID_SCORE;
	static const int INVALID_FRAME = -1;
	static const char INVALID_FLAG = '.';
	static const string INVALID_TOKEN;

	/* static methods */
	static VERSION guessVersion(const string& fn);
};

inline bool operator!=(const GFF& lhs, const GFF& rhs) {
	return !(lhs == rhs);
}

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* GFF_H_ */
