/*
 * SeqIO_test.cpp
 *
 *  Created on: Apr 25, 2018
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "PrimarySeq.h"
#include "SeqIO.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ostringstream;
using EGriceLab::MSGseqTK::PrimarySeq;
using EGriceLab::MSGseqTK::SeqIO;

int main() {
	PrimarySeq src1("ATCGNTCGANatcgntcgan", "seq1");
	PrimarySeq src2 = src1.revcom();

	/* try output */
	ostringstream out;
	SeqIO seqO(&out, "fasta");
	seqO.writeSeq(src1);
	seqO.setMaxLine(10); /* try multi-line output */
	seqO.writeSeq(src2);

	/* try input */
	string iStr = out.str();
	cout << "iStr: " << endl << iStr;

	istringstream in(iStr);
	SeqIO seqI(&in, "fasta");
	SeqIO stdO(&cout, "fasta");
	PrimarySeq dest1 = seqI.nextSeq();
	if(dest1 != src1)
		return 1;
	else
		stdO.writeSeq(dest1);

	PrimarySeq dest2 = seqI.nextSeq();
	if(dest2 != src2)
		return 1;
	else
		stdO.writeSeq(dest2);

	stdO.reset(&cout, "fastq");
	stdO.writeSeq(dest2);
}
