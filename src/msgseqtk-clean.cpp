/*******************************************************************************
 * This file is part of MSGseqTK, a Metagenomics Shot-Gun sequencing ToolKit
 * for ultra-fast and accurate MSG-seq cleaning, mapping and more,
 * based on space-efficient FM-index on entire collection of meta-genomics sequences.
 * Copyright (C) 2018  Qi Zheng
 *
 * MSGseqTK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MSGseqTK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
