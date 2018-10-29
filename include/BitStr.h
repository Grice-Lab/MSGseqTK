/*
 * BitStr.h
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#ifndef BITSTR_H_
#define BITSTR_H_

#include <string>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cassert>
#include "libsdsConst.h"
#include "libsdsBitBasic.h"

namespace EGriceLab {
namespace libSDS {

using std::basic_string;
using std::string;
using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

/**
 * A BitStr is a dynamic array representation of bits for arbitrary unsigned integer types
 */
template<typename uIntType = uint32_t>
class BitStr {
public:
	typedef uIntType value_type;
	typedef size_t size_type;
	/* constructors */
	/** default constructor */
	BitStr() = default;

	/**
	 * construct a BitStr with n-bits all set to zero
	 */
	explicit BitStr(size_type nB) : wid(sizeof(value_type) * Wb), nB(nB), n(int_len(nB, wid)), data(n, 0)
	{ 	}

	/**
	 * construct a BitStr with given length and a value,
	 * val will be used to fill the lowest/least-significant element of BitStr, with all others as zero
	 * @param nB  length in bits
	 * @param val  value to fill the lowest bits
	 */
	BitStr(size_type nB, uIntType val) : BitStr(nB) /* deligating construction */ {
		data[0] = val;
	}

	/**
	 * construct a BitStr by copying from a given C-type string,
	 * @param str  C-type string
	 * @param n  length of str
	 */
	BitStr(const uIntType* src, size_type n) : wid(sizeof(uIntType) * Wb), n(n), data(src, n) {
		nB = n * wid;
	}

	/**
	 * construct a BitStr by copying from a std__basic_string
	 * @param str  std__basic_string
	 */
	explicit BitStr(const basic_string<uIntType>& src) : wid(sizeof(uIntType) * Wb), nB(src.length() * wid), n(src.length()), data(src)
	{ 	}

	/**
	 * construct a BitStr by copying from a given C-type string of a diffent type,
	 * @param str  C-type string
	 * @param nSrc  length of str
	 */
	template<typename oIntType>
	BitStr(const oIntType* src, size_type nSrc) : wid(sizeof(uIntType) * Wb), nB(nSrc * wid),  n(int_len(nB, wid)), data(n, 0)
	{
		size_t wSrc = sizeof(oIntType) * Wb; // input width
		data.resize(n);
		if(wSrc < wid) {
			for(size_type i = 0; i < nB; i += wSrc)
				set(i, wSrc, src[i]);
		}
		else if(wSrc > wid) {
			for(size_type i = 0; i < nB; i += wid)
				setValue(i, src[i / wSrc] & (~0UL) << i % wSrc);
		}
		else
			std::copy(src, src + nSrc, data.begin());
	}

	/**
	 * construct a BitStr by copying from a std::basic_string a diffent type,
	 * @param str  std::basic_string
	 */
	template<typename oIntType>
	explicit BitStr(const basic_string<oIntType>& src) : wid(sizeof(uIntType) * Wb), nB(src.length() * wid), n(int_len(nB, wid)), data(n, 0)
	{
		size_t wSrc = sizeof(oIntType) * Wb; // input width
		if(wSrc < wid) {
			for(size_type i = 0; i < nB; i += wSrc)
				set(i, wSrc, src[i]);
		}
		else if(wSrc > wid) {
			for(size_type i = 0; i < nB; i += wid)
				setValue(i, src[i / wSrc] & (~0UL) << i % wSrc);
		}
		else
			std::copy(src.begin(), src.end(), data.begin());
	}

	/**
	 * construct a BitStr by coping from another BitStr of different type
	 * this non-explicit constructor can be used during copy-initiation as a converting constructor
	 */
	template<typename oIntType>
	BitStr(const BitStr<oIntType>& other) : wid(sizeof(value_type) * Wb), nB(other.length()), n(int_len(nB, wid)) {
		const size_t wO = other.getWid();
		data.resize(n);
		/* block copy */
		if(wO < wid) {
			for(size_type i = 0; i < nB; i += wO)
				set(i, wO, other.get(i, wO));
		}
		else if(wO > wid) {
			for(size_type i = 0; i < nB; i += wid)
				setValue(i, other.get(i, wid));
		}
		else {
			const basic_string<oIntType>& oData = other.getData();
			std::copy(oData.begin(), oData.end(), data.begin());
		}
	}

	/* member methods */
	/** get wid */
	uint getWid() const {
		return wid;
	}

	/** get length in bits */
	size_t length() const {
		return nB;
	}

	/** get number of bits, alias of length() */
	size_type numBits() const {
		return length();
	}

	/** get number of values */
	size_type numValues() const {
		return n;
	}

	/** test whether this BitStr is empty */
	bool empty() const {
		return n == 0;
	}

	/** resize this BitStr to nB bits, if the nB > length(), the remaining values are set to val,
	 * otherwise they are discarded
	 */
	void resize(size_type nB, bool val = false) {
		assert(wid > 0);
		if(nB == this->nB)
			return;
		size_type n = int_len(nB, wid);
		data.resize(n);

		/* fix last/lowest block of remaining bits, if any */
		for(size_t i = nB; i < n * wid; ++i)
			set(i, val);

		this->nB = nB;
		this->n = n;
	}

	/** clear all bits to zero */
	void clear() {
		std::fill(data.begin(), data.end(), 0);
	}

	/** reset this BitStr to default state */
	void reset() {
		nB = 0;
		n = 0;
		data.clear();
	}

	/** get the underlying data in raw array */
	const basic_string<uIntType>& getData() const {
		return data;
	}

	/**
	 * get the p-th element measured by value_type
	 */
	value_type getValue(size_type p) const {
		return data[p];
	}

	/**
	 * get the value represented by a given fixed length of bits
	 * @param start  start position measured by len
	 * @param len  number of bits
	 * @return value represented by these bits
	 */
	size_t getValue(size_type start, size_type len) const {
		assert(len <= wid);
		if(len == 0)
			return 0;
		if((start + 1) * len > nB)
			len = nB - start * len;
		size_type i = start * len / wid;
		size_type j = start * len - wid * i;
		size_t result;
		if (j + len <= wid)
			result = (data[i] << wid - j - len) >> (wid - len);
		else {
			result = data[i] >> j;
			result |= (data[i + 1] << wid * 2 - j - len) >> (wid - len);
		}
		return result;
	}

	/**
	 * set the n-th element to v measured by value_type
	 */
	void setValue(size_type n, value_type v) {
		data[n] = v;
	}

	/**
	 * set a given fixed region of fixed length bits to given value
	 * @param start  start position measured by size_type
	 * @param len  number of bits
	 * @param v  value to set
	 */
	void setValue(size_type start, size_type len, value_type v) {
		assert(len <= wid);
		if(len == 0)
			return;
		if((start + 1) * len > nB)
			len = nB - start * len;
		size_type i = start * len / wid;
		size_type j = start * len - i * wid;
		size_t mask = ((j + len) < wid ? ~0UL << j + len : 0UL)
			| ((wid - j) < wid ? ~0UL >> wid - j : 0UL);
		data[i] = (data[i] & mask) | v << j;
		if (j + len > wid) {
			mask = (~0UL) << len + j - wid;
			data[i+1] = (data[i + 1] & mask) | v >> wid - j;
		}
	}

	/* bit-wise methods */
	/** get the n-th bit of this BitStr */
	bool get(size_type i) const {
		return data[i / wid] & (1UL << (i % wid));
	}

	/**
	 * get the value represented by a region [start, start + len) of bits
	 * @param start  start position in bits
	 * @param len  region length in bits
	 */
	size_t get(size_type start, size_type len) const {
		assert(len <= wid);
		if(len == 0)
			return 0;
		if(start + len > nB)
			len = nB - start;
		size_t i = start / wid;
		size_t j = start - wid * i;
		size_t result;
		if (j + len <= W)
			result = (data[i] << (wid - j - len)) >> (wid - len);
		else {
			result = data[i] >> j;
			result |= (data[i+1] << (2 * wid - j - len)) >> (wid - len);
		}
		return result;
	}

	/** set the i-th bit of this BitStr */
	BitStr<uIntType>& set(size_type i, bool bit = true) {
		/* clear bits first */
		data[i / wid] &= ~(1UL << (i % wid));
		/* set bit */
		data[i / wid] |= bit << (i % wid);
		return *this;
	}

	/**
	 * set a given region of bits to a value
	 * @param start  start position in bits
	 * @param len  length in bits
	 * @param val  value to be set
	 */
	BitStr<uIntType>& set(size_type start, size_type len, size_t value) {
		assert(len <= wid);
		if(len == 0)
			return *this;
		if(start + len > nB)
			len = nB - start;
		size_t i = start / wid;
		size_t j = start - i * wid;
		size_t mask = ((j + len) < wid ? ~0UL << j + len : 0UL)
			| ((wid - j) < wid ? ~0UL >> wid - j : 0UL);
		data[i] = (data[i] & mask) | value << j;
		if (j + len > wid) {
			mask = (~0UL) << len + j - wid;
			data[i+1] = (data[i + 1] & mask) | value >> wid - j;
		}
		return *this;
	}

	/** reset the i-th bit to zero, alias of set(n, false) */
	BitStr<uIntType>& reset(size_type i) {
		data[i / wid] &= ~(1UL << (i % wid));
		return *this;
	}

	/** flip/toggle the i-th bit */
	BitStr<uIntType>& flip(size_type i) {
		data[i / wid] ^= 1UL << (i % wid);
		return *this;
	}

	/** test whether the i-th bit is on, alias to get(n) */
	bool test(size_type i) const {
		return get(i);
	}

	/** test whether any bit is on */
	bool any() const {
		return std::any_of(data.begin(), data.end(), [](uIntType v) { return v; });
	}

	/** test whether none of the bit is on */
	bool none() const {
		return !any();
	}

	/** test whether all bits are on */
	bool all() const {
		if(std::any_of(data.begin(), data.end() - 1, [](uIntType v) { return ~v; }))
			return false;
		for(size_type i = (n - 1) * wid; i < nB; ++i) /* test last block bits */
			if(!test(i))
				return false;
		return true;
	}

	/** count on bits using fast popcount */
	size_type count() const {
		if(empty())
			return 0;
		size_type on = 0;
		for(size_t b = 0; b < n - 1; ++b) /* popcount till last block */
			on += popcount(getValue(b));
		/* count the lowest bits bit by bit */
		for(size_t i = (n - 1) * wid; i < nB; ++i)
			if(test(i))
				on++;
		return on;
	}

	/**
	 * get a binary string representation
	 * @return  a binarhy string copy of this BitStr. A character in string is '1' if set, and '0' if not.
	 * Character position i in string corresponds to bit position b.size() - 1 - 1, in other word, reversed order
	 */
	string to_string() const {
		string bin;
		bin.reserve(length());
		for(size_type i = 0; i < nB; ++i)
			bin.push_back(test(i) ? '1' : '0');
		std::reverse(bin.begin(), bin.end());
		return bin;
	}

	/** write this BitStr to text output */
	ostream& write(ostream& out) const {
		out << to_string();
		return out;
	}

	/** save BitStr to binary output */
	ostream& save(ostream& out) const {
		out.write((const char*) &wid, sizeof(size_t));
		out.write((const char*) &nB, sizeof(size_type));
		out.write((const char*) &n, sizeof(size_type));
		out.write((const char*) data.c_str(), sizeof(uIntType) * n);
		return out;
	}

	/** load BitStr from binary input */
	istream& load(istream& in) {
		in.read((char*) &wid, sizeof(size_t));
		in.read((char*) &nB, sizeof(size_type));
		in.read((char*) &n, sizeof(size_type));
		uIntType *tmp = new uIntType[n];
		in.read((char*) tmp, sizeof(uIntType) * n);
		data.assign(tmp, n);
		delete[] tmp;
		return in;
	}

	/** test whether two BitStr are equal */
	template<typename oIntType>
	friend bool operator==(const BitStr<oIntType>& lhs, const BitStr<oIntType>& rhs);

	/** get required storage size of BitStr in bytes */
	size_t getBytes() const {
		return sizeof(wid) + sizeof(nB) + sizeof(n), sizeof(data) + sizeof(this);
	}

	/* member fields */
private:
	size_t wid = sizeof(uIntType) * Wb;    /* bit-width of a value_type */
	size_type nB = 0; /* number of bits */
	size_type n = 0; /* number of values in value_type */
	basic_string<uIntType> data; /* underlying data stored in std::basic_string, with length n */
};

typedef BitStr<uint8_t> BitStr8;
typedef BitStr<uint16_t> BitStr16;
typedef BitStr<> BitStr32;
typedef BitStr<uint64_t> BitStr64;

template<typename T>
inline bool operator==(const BitStr<T>& lhs, const BitStr<T>& rhs) {
	return lhs.wid == rhs.wid && lhs.nB == rhs.nB &&
			lhs.n == rhs.n && lhs.data == rhs.data;
}

/** test whether two BitStr are different */
template<typename T>
inline bool operator!=(const BitStr<T>& lhs, const BitStr<T>& rhs) {
	return !(lhs == rhs);
}

/** formatted write this BitStr as binary string */
template<typename T>
ostream& operator<<(ostream& out, const BitStr<T>& bstr) {
	bstr.write(out);
	return out;
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSTR_H_ */
