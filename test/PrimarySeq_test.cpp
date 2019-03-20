/*
 * PrimarySeq_test.cpp
 *
 *  Created on: Apr 25, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <algorithm>
#include "PrimarySeq.h"

using namespace std;
using namespace EGriceLab::MSGseqTK;

int main() {
	PrimarySeq src1("ATCGNTCGANatcgntcgan", "seq1");
	cout << "src1:" << endl << src1.getSeq() << endl;
	PrimarySeq dest1 = src1;

	if(dest1 != PrimarySeq("ATCGNTCGANatcgntcgan", "seq1"))
		return 1;
	else
		cout << "src1:" << endl << dest1.getSeq() << endl;

	dest1 = src1.revcom();
	if(dest1 != PrimarySeq("ntcgancgatNTCGANCGAT", "seq1"))
		return 1;
	else
		cout << "src1.revcom():" << endl << dest1.getSeq() << endl;
}
