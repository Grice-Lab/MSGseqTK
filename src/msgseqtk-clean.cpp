//============================================================================
// Name        : msgseqtk.cpp
// Author      : Qi Zheng
// Version     :
// Copyright   : GPL v3.0 Copyright (C) 2017  Qi Zheng
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "MetaGenomics Shotgun Sequencing cleaning by removing host contamination reads,"
		 << " based on Reduced-FM-index (RFM-index) powered Maximal Exact Matched Seeds (MEMS) searches"
		 << " and Baysian inference of background/target origin using multinomial models" << endl;
}

int main() {
	cout << "Hello World" << endl; // prints Hello World
	return 0;
}
