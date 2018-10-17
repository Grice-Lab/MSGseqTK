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

const BitSeqRRR::TableOffset BitSeqRRR::OFFSET(BitSeqRRR::BLOCK_SIZE); /* pre-computed TalbeOffset given block size */

BitSeqRRR::TableOffset::~TableOffset() {
	delete[] bitmaps;
	delete[] offset_class;
	delete[] rev_offset;
	for(uint32_t i = 0; i <= u; ++i) {
		delete[] binomial[i];
		delete[] log2binomial[i];
	}
	delete[] binomial;
	delete[] log2binomial;
}

BitSeqRRR::TableOffset::TableOffset(size_t u) : u(u) {
	bitmaps = new uint32_t[numBitmaps() + 1];
	rev_offset = new uint32_t[numRevOffset()];
	offset_class = new uint32_t[numClasses() + 1];
	binomial = new uint32_t*[u + 1];
	log2binomial = new uint32_t*[u + 1];
	for(uint32_t i = 0; i <= u; ++i) {
		binomial[i] = new uint32_t[u + 1](); /* value initiation */
		log2binomial[i] = new uint32_t[u + 1](); /* value initiation */
	}

	init_binomials();
	init_offsets();
}

BitSeqRRR::TableOffset::TableOffset(const TableOffset& other) : u(other.u) {
	bitmaps = new uint32_t[numBitmaps() + 1];
	rev_offset = new uint32_t[numRevOffset()];
	offset_class = new uint32_t[numClasses() + 1];
	binomial = new uint32_t*[u + 1];
	log2binomial = new uint32_t*[u + 1];
	for(uint32_t i = 0; i <= u; ++i) {
		binomial[i] = new uint32_t[u + 1](); /* value initiation */
		log2binomial[i] = new uint32_t[u + 1](); /* value initiation */
	}
	/* copy values */
	std::copy(other.bitmaps, other.bitmaps + other.numBitmaps() + 1, bitmaps);
	std::copy(other.rev_offset, other.rev_offset + other.numRevOffset(), rev_offset);
	std::copy(other.offset_class, other.offset_class + other.numClasses() + 1, offset_class);
	for(uint32_t i = 0; i <= u; ++i) {
		std::copy(other.binomial[i], other.binomial[i] + u + 1, binomial[i]);
		std::copy(other.log2binomial[i], other.log2binomial[i] + u + 1, log2binomial[i]);
	}
}

void BitSeqRRR::TableOffset::init_binomials() {
	for(uint32_t i = 0; i <= u; ++i) { /* upper-right half all ones */
		binomial[i][0] = 1;
		binomial[i][1] = 1;
		binomial[i][i] = 1;
		log2binomial[i][0] = 0;
		log2binomial[i][1] = 0;
		log2binomial[i][i] = 0;
	}
	for(uint32_t j = 1; j <= u; ++j) {
		for(uint32_t i = j + 1; i <= u; ++i) {
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
	for (uint32_t i = start; i < u; ++i)
		idx += init_classes(shift, classIdx, k, len + 1, i + 1, val | (1UL << i)); /* recursive initiation */
	return idx;
}

void BitSeqRRR::TableOffset::init_offsets() {
	uint32_t shift = 0;
	uint32_t classIdx = 0;
	offset_class[0] = 0;
	for (uint32_t k = 0; k <= u; ++k) {
		shift += init_classes(shift, classIdx, k);
		offset_class[k + 1] = classIdx;
	}
}

void BitSeqRRR::build_sampled() {
	/* sampling C */
	nCsampled = numClassSampled();
	wCsampled = bits(ones);
	Csampled = BitStr<uint32_t>(nCsampled * wCsampled);
	size_t sum = 0;
	for(size_t i = 0; i < nC; ++i) {
		if(i % sample_rate == 0)
			Csampled.setValue(i / sample_rate, wCsampled, sum);
		sum += C.getValue(i, wC);
	}
	/* set the last field */
	for(size_t i = (nC + sample_rate - 1) / sample_rate; i < nCsampled; ++i) {
		Csampled.setValue(i, wCsampled, sum);
	}

	/* sampling O */
	nOsampled = numOffsetSampled();
	wOsampled = bits(nO);
	Osampled = BitStr<uint32_t>(nOsampled * wOsampled);
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
	while(k < pos - 1 && pos > 0) {
		size_t lower = C.getValue(k, wC) & mask;
		size_t upper = C.getValue(k + 1, wC);
		sum += lower + upper;
		posO += OFFSET.get_log2binomial(BLOCK_SIZE, lower) + OFFSET.get_log2binomial(BLOCK_SIZE, upper);
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
	while(start < nC - 1 && acc == Csampled.getValue(start + 1, wCsampled))
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
	assert(acc == r);
	assert(rank1(pos) == r);
	assert(access(pos));
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
	while(start < end - 1) {
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
	while(start < nC -1 && acc + sample_rate * BLOCK_SIZE == Csampled.getValue(start + 1, wCsampled)) {
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
	assert(acc == r);
	assert(rank0(pos) == r);
	assert(!access(pos));
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

BitSeqRRR& BitSeqRRR::swap(BitSeqRRR& other) {
	C.swap(other.C);
	O.swap(other.O);
	Csampled.swap(other.Csampled);
	Osampled.swap(other.Osampled);

	std::swap(nC, other.nC);
	std::swap(wC, other.wC);
	std::swap(nO, other.nO);
	std::swap(wO, other.wO);

	std::swap(nCsampled, other.nCsampled);
	std::swap(wCsampled, other.wCsampled);
	std::swap(nOsampled, other.nOsampled);
	std::swap(wOsampled, other.wOsampled);

	std::swap(sample_rate, other.sample_rate);
	return *this;
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
	/* compare basic fields */
	if(dynamic_cast<const BitSeq&>(lhs) != dynamic_cast<const BitSeq&>(rhs))
		return false;
	cerr << "basic equals" << endl;
	/* compare member fields */
	if(lhs.sample_rate != rhs.sample_rate)
		return false;
	cerr << "rate equals" << endl;
	if(!(lhs.nC == rhs.nC && lhs.wC == rhs.wC && lhs.C == rhs.C))
		return false;
	cerr << "C equals" << endl;
	if(!(lhs.nO == rhs.nO && lhs.wO == rhs.wO && lhs.O == rhs.O))
		return false;
	cerr << "O equals" << endl;

	if(!(lhs.nCsampled == rhs.nCsampled && lhs.wCsampled == rhs.wCsampled && lhs.Csampled == rhs.Csampled))
		return false;
	cerr << "Csampled equals" << endl;

	if(!(lhs.nOsampled == rhs.nOsampled && lhs.wOsampled == rhs.wOsampled && lhs.Osampled == rhs.Osampled))
		return false;
	cerr << "Osampled equals" << endl;

	cerr << "returning true" << endl;
	return true;
}

} /* namespace libSDS */
} /* namespace EGriceLab */
