/*
 * BitSeqRRR.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: zhengqi
 */

#include <algorithm>
#include "BitSeqRRR.h"

namespace EGriceLab {
namespace libSDS {

const BitSeqRRR::TableOffset BitSeqRRR::OFFSET; /* pre-computed TalbeOffset with pre-defined BLOCK_SIZE */

BitSeqRRR::TableOffset::TableOffset() {
	init_binomials();
	init_offsets();
}

void BitSeqRRR::TableOffset::init_binomials() {
	for(uint32_t i = 0; i <= BLOCK_SIZE; ++i) { /* upper-right half all ones */
		binomial[i][0] = 1;
		binomial[i][1] = 1;
		binomial[i][i] = 1;
		log2binomial[i][0] = 0;
		log2binomial[i][1] = 0;
		log2binomial[i][i] = 0;
	}
	for(uint32_t j = 1; j <= BLOCK_SIZE; ++j) {
		for(uint32_t i = j + 1; i <= BLOCK_SIZE; ++i) {
			binomial[i][j] = binomial[i - 1][j - 1] + binomial[i - 1][j];
			log2binomial[i][j] = bits(binomial[i][j] - 1);
		}
	}
}

uint32_t BitSeqRRR::TableOffset::init_classes(uint32_t& shift, uint32_t& classIdx, uint32_t k,
		uint32_t len, uint32_t start, uint32_t val) {
	uint idx = 0;
	if (k == len) {
		bitmaps[classIdx] = val;
		rev_offset[val] = classIdx - shift;
		classIdx++;
		return 1;
	}
	if (k < len)
		return 0;
	for (uint32_t i = start; i < BLOCK_SIZE; ++i)
		idx += init_classes(shift, classIdx, k, len + 1, i + 1, val | (1UL << i)); /* recursive initiation */
	return idx;
}

void BitSeqRRR::TableOffset::init_offsets() {
	uint32_t shift = 0;
	uint32_t classIdx = 0;
	offset_class[0] = 0;
	for (uint32_t k = 0; k <= BLOCK_SIZE; ++k) {
		shift += init_classes(shift, classIdx, k);
		offset_class[k + 1] = classIdx;
	}
}

void BitSeqRRR::build_sampled() {
	/* sampling C */
	nCsampled = numClassSampled();
	wCsampled = bits(ones);
	Csampled.resize(nCsampled * wCsampled);
	size_t sum = 0;
	for(size_t i = 0; i < nC; ++i) {
		if(i % sample_rate == 0)
			Csampled.setValue(i / sample_rate, wCsampled, sum);
		sum += C.getValue(i, wC);
	}
	/* set the last field */
	for(size_t i = (nC + sample_rate - 1) / sample_rate; i < nCsampled; ++i)
		Csampled.setValue(i, wCsampled, sum);

	/* sampling O */
	nOsampled = numOffsetSampled();
	wOsampled = bits(nO);
	Osampled.resize(nOsampled * wOsampled);
	for(size_t i = 0, pos = 0; i < nC; ++i) {
		if(i % sample_rate == 0)
			Osampled.setValue(i / sample_rate, wOsampled, pos);
		pos += OFFSET.get_log2binomial(BLOCK_SIZE, C.getValue(i, wC));
	}
}

size_t BitSeqRRR::getBytes() const {
	return BitSeq::getBytes() +
			C.getBytes() + O.getBytes() + sizeof(nC) + sizeof(nO) + sizeof(wC) + sizeof(wO) +
			Csampled.getBytes() + sizeof(nCsampled) + sizeof(wCsampled) +
			Osampled.getBytes() + sizeof(nOsampled) + sizeof(wOsampled) +
			sizeof(sample_rate) + sizeof(this);
}

size_t BitSeqRRR::rank1(size_t i) const {
	size_t nearest_sampled_value = i / BLOCK_SIZE / sample_rate;
	size_t sum = Csampled.getValue(nearest_sampled_value, wCsampled);
	size_t posO = Osampled.getValue(nearest_sampled_value, wOsampled);
	size_t pos = i / BLOCK_SIZE;
	size_t k = nearest_sampled_value * sample_rate;
	if(k % 2 == 1 && k < pos) {
		size_t aux = C.getValue(k, wC);
		sum += aux;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, aux);
		k++;
	}
	size_t mask = 0x0F;
	const uint8_t* arr = reinterpret_cast<const uint8_t*>(C.getData().c_str());
	arr += k/2;
	while(k + 1 < pos) {
//		size_t lower = C.getValue(k, wC) & mask;
//		size_t upper = C.getValue(k + 1, wC);
		size_t lower = *arr & mask;
		size_t upper = *arr / 16;
//		assert(lower == C.getValue(k, wC) & mask);
//		assert(upper == C.getValue(k + 1, wC));
		sum += lower + upper;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, lower) + OFFSET.get_log2binomial(BLOCK_SIZE, upper);
		arr++;
		k += 2;
	}
	if(k < pos) { /* process last field */
		size_t aux = C.getValue(k, wC);
		sum += aux;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, aux);
		k++;
	}
	size_t c = C.getValue(pos, wC);
	sum += popcount32(((2UL << (i % BLOCK_SIZE)) - 1) &
			OFFSET.get_bitmap(c, O.get(posO, OFFSET.get_log2binomial(BLOCK_SIZE, c))));
	return sum;
}

size_t BitSeqRRR::select1(size_t r) const {
	if(r == 0)
		return -1;
	if(r > ones)
		return n;
	// Search over partial sums
	size_t start = 0;
	size_t end = nCsampled - 1;
	size_t med, acc = 0, pos;
	while(start < end - 1) {
		med = (start + end) / 2;
		acc = Csampled.getValue(med, wCsampled);
		if(acc < r) {
			if(med == start)
				break;
			start = med;
		}
		else {
			if(end == 0)
				break;
			end = med-1;
		}
	}
	acc = Csampled.getValue(start, wCsampled);
	while(start + 1 < nC && acc == Csampled.getValue(start + 1, wCsampled))
		start++;
	pos = start * sample_rate;
	size_t posO = Osampled.getValue(start, wOsampled);
	acc = Csampled.getValue(start, wCsampled);

	// Sequential search over C
	size_t s;
	for(s = 0; pos < nC; ++pos) {
		s = C.getValue(pos, wC);
		if(acc + s >= r)
			break;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, s);
		acc += s;
	}
	pos *= BLOCK_SIZE;

	// Search inside the block
	while(acc < r) {
		size_t new_posO = posO + OFFSET.get_log2binomial(BLOCK_SIZE, s);
		size_t block = OFFSET.get_bitmap(s, O.get(posO, new_posO - posO));
		posO = new_posO;
		new_posO = 0;
		while(acc < r && new_posO < BLOCK_SIZE) {
			pos++;
			new_posO++;
			acc += (block & 1) != 0;
			block /= 2;
		}
	}
	pos--;
//	assert(acc == r);
//	assert(rank1(pos) == r);
//	assert(access(pos));
	return pos;
}

size_t BitSeqRRR::select0(size_t r) const {
	if(r == 0)
		return -1;
	if(r > numZeros())
		return n;

	// Search over partial sums
	size_t start = 0;
	size_t end = nCsampled - 1;
	size_t med, acc = 0, pos;
	while(start + 1 < end) {
		med = (start + end) / 2;
		acc = med * sample_rate * BLOCK_SIZE - Csampled.getValue(med, wCsampled);
		if(acc < r) {
			if(med == start)
				break;
			start = med;
		}
		else {
			if(end == 0)
				break;
			end = med-1;
		}
	}
	acc = Csampled.getValue(start, wCsampled);
	while(start + 1 < nC && acc + sample_rate * BLOCK_SIZE == Csampled.getValue(start + 1, wCsampled)) {
		start++;
		acc += sample_rate * BLOCK_SIZE;
	}
	acc = start * sample_rate * BLOCK_SIZE - acc;
	pos = start * sample_rate;
	size_t posO = Osampled.getValue(start, wOsampled);

	// Sequential search over C
	size_t s;
	for(s = 0; pos < nC; ++pos) {
		s = C.getValue(pos, wC);
		if(acc + BLOCK_SIZE - s >=r)
			break;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, s);
		acc += BLOCK_SIZE - s;
	}
	pos *= BLOCK_SIZE;

	// Search inside the block
	while(acc < r) {
		size_t new_posO = posO + OFFSET.get_log2binomial(BLOCK_SIZE, s);
		size_t block = OFFSET.get_bitmap(s, O.get(posO, new_posO - posO));
		posO = new_posO;
		new_posO = 0;
		while(acc < r && new_posO < BLOCK_SIZE) {
			pos++;
			new_posO++;
			acc += (block & 1) == 0;
			block /= 2;
		}
	}
	pos--;
//	assert(acc == r);
//	assert(rank0(pos) == r);
//	assert(!access(pos));
	return pos;
}

bool BitSeqRRR::access(size_t i) const {
	size_t nearest_sampled_value = i / BLOCK_SIZE / sample_rate;
	size_t posO = Osampled.getValue(nearest_sampled_value, wOsampled);
	size_t pos = i / BLOCK_SIZE;
	assert(pos <= nC);
	for(size_t k = nearest_sampled_value * sample_rate; k < pos; ++k)
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, C.getValue(k, wC));
	size_t c = C.getValue(pos, wC);
	return ((1UL << (i % BLOCK_SIZE)) &
			OFFSET.get_bitmap(c, O.get(posO, OFFSET.get_log2binomial(BLOCK_SIZE, c)))) != 0;
}

bool BitSeqRRR::access(size_t i, size_t& r) const {
	if(i == -1)
		return 0;
	size_t nearest_sampled_value = i / BLOCK_SIZE / sample_rate;
	size_t sum = Csampled.getValue(nearest_sampled_value, wCsampled);
	size_t posO = Osampled.getValue(nearest_sampled_value, wOsampled);
	size_t pos = i / BLOCK_SIZE;
	size_t k = nearest_sampled_value * sample_rate;
	if(k % 2 == 1 && k < pos) {
		size_t aux = C.getValue(k, wC);
		sum += aux;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, aux);
		k++;
	}
	size_t mask = 0x0F;
	const uint8_t* arr = reinterpret_cast<const uint8_t*>(C.getData().c_str());
	arr += k / 2;
	while(k + 1 < pos) {
//		size_t lower = C.getValue(k, wC) & mask;
//		size_t upper = C.getValue(k + 1, wC);
		size_t lower = *arr & mask;
		size_t upper = *arr / 16;
		sum += lower + upper;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, lower) + OFFSET.get_log2binomial(BLOCK_SIZE, upper);
		arr++;
		k += 2;
	}
	if(k < pos) {
		size_t aux = C.getValue(k, wC);
		sum += aux;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, aux);
		k++;
	}
	size_t c = C.getValue(pos, wC);
	size_t v = OFFSET.get_bitmap(c, O.get(posO, OFFSET.get_log2binomial(BLOCK_SIZE, c)));
	sum += popcount32(((2UL << (i % BLOCK_SIZE)) - 1) & v);
	r = sum;
	if( ((1UL << (i % BLOCK_SIZE)) & v))
		return true;
	else {
		r = i - r + 1;
		return false;
	}
}

ostream& BitSeqRRR::save(ostream& out) const {
	BitSeq::save(out); /* save base object */
	assert(C.length() == nC * wC);
	assert(O.length() == nO * wO);
	assert(Csampled.length() == nCsampled * wCsampled);
	assert(Osampled.length() == nOsampled * wOsampled);

	/* save C and O */
	C.save(out);
	O.save(out);
	out.write((const char*) &nC, sizeof(size_t));
	out.write((const char*) &wC, sizeof(size_t));
	out.write((const char*) &nO, sizeof(size_t));
	out.write((const char*) &wO, sizeof(size_t));

	/* save sampled */
	Csampled.save(out);
	Osampled.save(out);
	out.write((const char*) &nCsampled, sizeof(size_t));
	out.write((const char*) &wCsampled, sizeof(size_t));
	out.write((const char*) &nOsampled, sizeof(size_t));
	out.write((const char*) &wOsampled, sizeof(size_t));

	out.write((const char*) &sample_rate, sizeof(size_t));
	return out;
}

istream& BitSeqRRR::load(istream& in) {
	BitSeq::load(in); /* load base object */
	/* load C and O */
	C.load(in);
	O.load(in);
	in.read((char*) &nC, sizeof(size_t));
	in.read((char*) &wC, sizeof(size_t));
	in.read((char*) &nO, sizeof(size_t));
	in.read((char*) &wO, sizeof(size_t));

	/* load sampled */
	Csampled.load(in);
	Osampled.load(in);
	in.read((char*) &nCsampled, sizeof(size_t));
	in.read((char*) &wCsampled, sizeof(size_t));
	in.read((char*) &nOsampled, sizeof(size_t));
	in.read((char*) &wOsampled, sizeof(size_t));

	in.read((char*) &sample_rate, sizeof(size_t));

	assert(C.length() == nC * wC);
	assert(O.length() == nO * wO);
	assert(Csampled.length() == nCsampled * wCsampled);
	assert(Osampled.length() == nOsampled * wOsampled);
	return in;
}

bool operator==(const BitSeqRRR& lhs, const BitSeqRRR& rhs) {
	return dynamic_cast<const BitSeq&>(lhs) == dynamic_cast<const BitSeq&>(rhs) &&
			lhs.nC == rhs.nC && lhs.wC == rhs.wC && lhs.C == rhs.C &&
			lhs.nO == rhs.nO && lhs.wO == rhs.wO && lhs.O == rhs.O &&
			lhs.nCsampled == rhs.nCsampled && lhs.wCsampled == rhs.wCsampled && lhs.Csampled == rhs.Csampled &&
			lhs.nOsampled == rhs.nOsampled && lhs.wOsampled == rhs.wOsampled && lhs.Osampled == rhs.Osampled &&
			lhs.sample_rate == rhs.sample_rate;
}

} /* namespace libSDS */
} /* namespace EGriceLab */
