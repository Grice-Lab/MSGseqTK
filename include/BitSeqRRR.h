/*
 * BitSeqRRR.h
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#ifndef BITSEQRRR_H_
#define BITSEQRRR_H_

#include <utility>
#include "libsdsBitBasic.h"
#include "BitSeq.h"
#include "BitStr.h"

/*
 * block size can't be changed in this implementation
 * it would require more than just changing the constant
 */
#define BLOCK_SIZE (15)

namespace EGriceLab {
namespace libSDS {

/**
 * A Raman, Raman and Rao's implementation of BitSeq, based on the algorithm described in
 *  [1] R. Raman, V. Raman and S. Rao. Succinct indexable dictionaries with applications
 *     to encoding $k$-ary trees and multisets. SODA02.
 *  [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 */
class BitSeqRRR: public BitSeq {
public:
	/**
	 * pre-computed table offsets used exclusively for BitSequenceRRR algorithm
	 */
	class TableOffset {
	public:
		/* constructors */
		/** default constructor */
		TableOffset();

		/* member methods */
		/** computes binomial(n,k) for n,k <= u */
		uint32_t get_binomial(uint32_t n, uint32_t k) const {
			return binomial[n][k];
		}

		/** computes ceil(log2(binomial(n,k))) for n, k <= u */
		uint32_t get_log2binomial(uint32_t n, uint32_t k) const {
			return log2binomial[n][k];
		}

		/** get the bitmap represented by the given class and inclass offsets */
		uint32_t get_bitmap(uint32_t class_offset, uint inclass_offset) const {
			if(class_offset == 0)
				return 0;
			if(class_offset == BLOCK_SIZE)
				return (1 << BLOCK_SIZE) - 1;
			return bitmaps[offset_class[class_offset] + inclass_offset];
		}

		/** computes the offset of the first u bits of a given bitstring */
		uint32_t compute_offset(uint32_t v) const {
			return rev_offset[v];
		}

		/* utility methods */
	private:
		/** initiate binomial matrices */
		void init_binomials();

		/** initiate offset arrays */
		void init_offsets();

		/**
		 * initiate class arrays
		 * @param shift  class shift
		 * @param outer classIdx  class index
		 * @param k  inner class index
		 * @param len  max len
		 * @param start  start position
		 * @param val  initial value to use
		 */
		uint32_t init_classes(uint32_t& shift, uint32_t& classIdx, uint32_t k,
				uint32_t len = 0, uint32_t start = 0, uint32_t val = 0);

		/* static fields */
	public:
		static const uint32_t NUM_BITMAPS = 1 << BLOCK_SIZE;
		static const uint32_t NUM_REVOFFSETS = 2 << (BLOCK_SIZE + 1);
		static const uint32_t NUM_CLASSES = BLOCK_SIZE + 1;

		/* member fields */
	private:
		uint32_t binomial[BLOCK_SIZE + 1][BLOCK_SIZE + 1]; /* pre-computed binomial static matrix */
		uint32_t log2binomial[BLOCK_SIZE + 1][BLOCK_SIZE + 1]; /* pre-computed log2-binomial static matrix */
		uint32_t rev_offset[NUM_REVOFFSETS]; /* pre-computed reverse offset static array */
		uint32_t offset_class[NUM_CLASSES + 1]; /* pre-computed integer class-offset static array */
		uint32_t bitmaps[NUM_BITMAPS + 1]; /* pre-computed bitmap static array */
	};

	/* constructors */
	/** default constructor */
	BitSeqRRR() = default;

	/** construct a BitSeqRRR from a given BitStr with any type */
	template<typename oIntType>
	explicit BitSeqRRR(const BitStr<oIntType>& bstr, size_t sample_rate = DEFAULT_SAMPLE_RATE) : sample_rate(sample_rate) {
		build_basic(bstr);
		build_sampled();
	}

	/* member methods */
	/** get number of classes */
	size_t numClasses() const {
		return (n + BLOCK_SIZE - 1) / BLOCK_SIZE; /* ceil(n / BLOCK_SIZE) */
	}

	/** get number of sampled classes */
	size_t numClassSampled() const {
		return nC / sample_rate + 2;
	}

	/** get number of sampled offsets */
	size_t numOffsetSampled() const {
		return (nC + sample_rate - 1) / sample_rate; /* ceil(nC / sample_rate) */
	}

	/**
	 * get the size of the structure in bytes
	 * @override  base class virtual method
	 */
	virtual size_t getBytes() const;

	/**
	 * get # of ones until position i (0-based, inclusive)
	 * @param i  position
	 * @override  base class virtual method
	 * @return  number of ones until position i,
	 */
	virtual size_t rank1(size_t i) const;

	/**
	 * get the position of the i-th one
	 * @param r  the order/rank of the 1
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > ones
	 * @override  base class virtual method
	 */
	virtual size_t select1(size_t r) const;

	/**
	 * get the position of the i-th zero
	 * @param r  the order/rank of the 0
	 * @return  position in this BitSeq, or -1 if i = 0, or len if i > zeros
	 * @override  base class virtual method
	 * first binary search over first level rank structure
	 * then sequential search using popcount over a int
	 * then sequential search using popcount over a char
	 * then sequential search bit by bit
	 */
	virtual size_t select0(size_t r) const;

	/**
	 * get the i-th bit of this BitSeq by testing the rank1 difference
	 * @param i  pos
	 * @return  true if it is i-th bit is significant (one)
	 * @override  base class virtual method
	 */
	virtual bool access(size_t i) const;

	/**
	 * get the i-th bit of this BitSeq by testing the rank1 difference, and return the corresponding rank the same time
	 * @param i  pos
	 * @param r  rank of i
	 * @return  true if it is i-th bit is significant (one)
	 * @override  base class virtual method
	 */
	virtual bool access(size_t i, size_t& r) const;

	/**
	 * save this BitSeq to binary output
	 * @param out  binary output
	 * @override  base class virtual method
	 */
	virtual ostream& save(ostream& out) const;

	/**
	 * load data from a binary input
	 * @param in  binary input
	 * @override  base class virtual method
	 */
	virtual istream& load(istream& in);

	/**
	 * reset this BitSeqRRR to default state
	 * @override  base class method
	 */
	virtual void reset() {
		sample_rate = 0;
		nOsampled = wOsampled = 0;
		nCsampled = wCsampled = 0;
		Osampled.reset();
		Csampled.reset();
		nO = wO = 0;
		nC = wC = 0;
		O.reset();
		C.reset();
		BitSeq::reset();
	}

	/* utility method */
private:
	/** build the basic C and O data from any type */
	template<typename oIntType>
	void build_basic(const BitStr<oIntType>& bstr) {
		n = bstr.length();
		ones = 0;
		/* build Table C */
		wC = bits(BLOCK_SIZE);
		nC = numClasses();
		C.resize(nC * wC);

		nO = 0; /* number of total bits required for O */
		wO = 1; /* Offset in 1 bit */
		for(size_t i = 0; i < nC; ++i) {
			uint32_t value = popcount32(bstr.get(i * BLOCK_SIZE, BLOCK_SIZE));
//			assert(value <= BLOCK_SIZE);
			C.setValue(i, wC, value);
			ones += value;
			nO += OFFSET.get_log2binomial(BLOCK_SIZE, value);
		}

		/* build Table O */
		O.resize(nO * wO);
		for(size_t i = 0, Opos = 0; i < nC; ++i) {
			uint32_t value = bstr.get(i * BLOCK_SIZE, BLOCK_SIZE);
			O.set(Opos, OFFSET.get_log2binomial(BLOCK_SIZE, popcount32(value)), OFFSET.compute_offset(value));
			Opos += OFFSET.get_log2binomial(BLOCK_SIZE, popcount32(value));
		}
	}

	/** build the sampled C and O */
	void build_sampled();

	/* member fields */
private:
	BitStr32 C; /* bitstring for classes */
	BitStr32 O; /* bitstring for offsets */
	size_t nC = 0, wC = 0; /* len and bit wid for C */
	size_t nO = 0, wO = 0; /* len and bit wid for O */
	BitStr32 Csampled; /* C samplings */
	BitStr32 Osampled; /* O samplings */
	size_t nCsampled = 0, wCsampled = 0;
	size_t nOsampled = 0, wOsampled = 0;
	/** Sample rate */
	size_t sample_rate = 0;

	/* non-member methods */
	/* relationship operators */
public:
	friend bool operator==(const BitSeqRRR& lhs, const BitSeqRRR& rhs);

	/* static members */
public:
	static const size_t DEFAULT_SAMPLE_RATE = 32;
	static const TableOffset OFFSET; /* pre-computed TalbeOffset given block size */
};

inline bool operator!=(const BitSeqRRR& lhs, const BitSeqRRR& rhs) {
	return !(lhs == rhs);
}

} /* namespace libSDS */
} /* namespace EGriceLab */

#endif /* BITSEQRRR_H_ */
