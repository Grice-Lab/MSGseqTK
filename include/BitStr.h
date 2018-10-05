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
#include "libsdsConst.h"

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
template<typename uIntType = uint>
class BitStr {
public:
	typedef uIntType value_type;
	typedef size_t size_type;
	/* constructors */
	/** default constructor */
	BitStr() = default;

	/** copy constructor */
	BitStr(const BitStr<uIntType>& other) : wid(other.wid), n(other.n) {
		data = new uIntType[n];
		std::copy(other.data, other.data + other.n, data);
	}

	/** copy assignment operator using copy-swap */
	BitStr& operator=(BitStr<uIntType> other) {
		swap(other);
		return *this;
	}

	/** destructor */
	~BitStr() {
		delete[] data;
	}

	/**
	 * construct a BitStr with n-bits all set to zero
	 */
	BitStr(size_type nB) : wid(sizeof(value_type) * Wb) {
		n = (nB + wid - 1) / wid; /* ceil(nB / wid) */
		data = new uIntType[n]; /* value initialization */
	}

	/**
	 * construct a BitStr with given length and a value,
	 * val will be used to fill the least-significant element of BitStr, with all others as zero
	 * @param n  length of bits
	 * @param val  value to use this BitStr
	 */
	BitStr(size_type n, uIntType val) : wid(sizeof(value_type) * Wb), n(n) {
		data = new uIntType[n](); /* value initialization */
		data[0] = val;
	}

	/**
	 * construct a BitStr by copying from a given C-type string,
	 * @param str  C-type string
	 * @param n  length of str
	 */
	BitStr(const uIntType* src, size_type n) : wid(sizeof(value_type) * Wb), n(n) {
		data = new uIntType[n];
		std::copy(src, src + n, data);
	}

	/** construct a BitStr by coping from another BitStr of different type */
	template<typename oIntType>
	BitStr(const BitStr<oIntType>& other) : wid(sizeof(value_type) * Wb) {
		n = other.length() / wid;
		data = new uIntType[n];
		/* bitwise copy */
		for(size_type i = 0; i < length(); ++i)
			set(i, other.get(i));
	}

	/** copy assign operator by copying from another BitStr of potential different type */
	template<typename oIntType>
	BitStr<uIntType>& operator=(const BitStr<oIntType>& other) {
		BitStr<uIntType> newBstr(other);
		swap(newBstr);
		return *this;
	}

	/* member methods */
	/** get wid */
	uint getWid() const {
		return wid;
	}

	/** get number of bits */
	size_t numBits() const {
		return n * wid;
	}

	/** get length in bits, alias of numBits */
	size_type length() const {
		return numBits();
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
		size_type n = (nB + wid - 1) / wid;
		if(n == this->n)
			return;
		uIntType* data_new = new uIntType[n];
		if(data_new == nullptr) /* failed, do not change the data */
			return;
		if(n > this->n) {
			std::copy(data, data + this->n, data_new);
			if(val)
				std::fill(data + this->n, data + n, ~0);
		}
		else
			std::copy(data, data + n, data_new);

		std::swap(data, data_new); /* swap the array */
		this->n = n;
	}

	/** clear all bits to zero */
	void clear() {
		std::fill(data, data + n, 0);
	}

	/** get the underlying data in raw array */
	const uIntType* getData() const {
		return data;
	}

	/**
	 * get the p-th element measured by value_type
	 */
	value_type getValue(size_type p) const {
		return data[p];
	}

	/**
	 * set the p-th element to v measured by value_type
	 */
	void setValue(size_type p, value_type v) {
		data[p] = v;
	}

	/* bit-wise methods */
	/** get the n-th bit of this BitStr */
	bool get(size_type n) const {
		return data[n / wid] & (static_cast<uIntType>(1) << (n % wid));
	}

	/** set the n-th bit of this BitStr */
	BitStr<uIntType>& set(size_type n, bool bit = true) {
		/* clear bits first */
		data[n / wid] &= ~(static_cast<uIntType>(1) << (n % wid));
		/* set bit */
		data[n / wid] |= static_cast<uIntType>(bit) << (n % wid);
		return *this;
	}

	/** reset the n-th bit to zero, alias of set(n, true) */
	BitStr<uIntType>& reset(size_type n) {
		data[n / wid] &= ~(static_cast<uIntType>(1) << (n % wid));
		return *this;
	}

	/** flip/toggle the n-th bit */
	BitStr<uIntType>& flip(size_type n) {
		data[n / wid] ^= static_cast<uIntType>(1) << n;
		return *this;
	}

	/** test whether the n-th bit is on, alias to get(n) */
	bool test(size_type n) const {
		return get(n);
	}

	/** test whether any bit is on */
	bool any() const {
		for(size_type i = 0; i < n; ++i)
			if(data[i])
				return true;
		return false;
	}

	/** test whether none of the bit is on */
	bool none() const {
		return !any();
	}

	/** test whether all bits are on */
	bool all() const {
		for(size_type i = 0; i < n; ++i)
			if(~ data[i] != 0)
				return false;
		return true;
	}

	/** count on bits */
	size_type count() const {
		size_type on = 0;
		for(size_type i = 0; i < length(); ++i)
			if(test(i))
				on++;
		return on;
	}

	/** get a binary string representation */
	string to_string(const string& sep = "") const {
		string bin;
		bin.reserve(numBits());
		for(size_type i = 0; i < n; ++i) {
			if(!sep.empty() && i != 0)
				bin.append(sep);
			for(size_type j = wid; j != 0; --j) /* start from most-significant bit */
				bin.push_back(test(i * wid + j - 1) ? '1' : '0');
		}
		return bin;
	}

	/** write this BitStr to text output */
	ostream& write(ostream& out) const {
		out << to_string();
		return out;
	}

	/** test whether two BitStr are equal */
	template<typename oIntType>
	friend bool operator==(const BitStr<oIntType>& lhs, const BitStr<oIntType>& rhs);

	/** swap this BitStr with another with same type, the underlying data is shallow-swapped */
	void swap(BitStr<uIntType>& other) {
		std::swap(wid, other.wid);
		std::swap(n, other.n);
		std::swap(data, other.data);
	}

	/* member fields */
private:
	uint wid = W;    /* bit-width of a value_type */
	size_type n = 0; /* number of elements of value_type */
	uIntType* data = nullptr; /* underlying data of given type */
};

template<typename T>
inline bool operator==(const BitStr<T>& lhs, const BitStr<T>& rhs) {
	if(!(lhs.wid == rhs.wid && lhs.n == rhs.n))
		return false;
	for(size_t i = 0; i < rhs.n; ++i)
		if(lhs.data[i] != rhs.data[i])
			return false;
	return true;
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
