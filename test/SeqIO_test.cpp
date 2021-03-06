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

using namespace std;
using namespace EGriceLab::MSGseqTK;

int main() {
	PrimarySeq src1("ATCGNTCGANatcgntcgan", "seq1");
	PrimarySeq src2 = src1.revcom();

	/* try output */
	ostringstream out;
	SeqIO seqO(&out, SeqIO::FASTA);
	seqO.writeSeq(src1);
	seqO.setMaxLine(10); /* try multi-line output */
	seqO.writeSeq(src2);

	/* try input */
	string iStr = out.str();
	cout << "iStr: " << endl << iStr;

	istringstream in(iStr);
	SeqIO seqI(&in, SeqIO::FASTA);
	SeqIO stdO(&cout, SeqIO::FASTA);
	PrimarySeq dest1 = seqI.nextSeq();
	if(dest1 != src1) {
		cerr << "Read in dest1 is different than src1" << endl <<
				"src1:  " << src1.getSeq() << endl <<
				"dest1: " << dest1.getSeq() << endl;
		return EXIT_FAILURE;
	}
	else
		stdO.writeSeq(dest1);

	PrimarySeq dest2 = seqI.nextSeq();
	if(dest2 != src2)
		return 1;
	else
		stdO.writeSeq(dest2);

	stdO.reset(&cout, SeqIO::FASTQ);
	stdO.writeSeq(dest2.fixQual());
}
