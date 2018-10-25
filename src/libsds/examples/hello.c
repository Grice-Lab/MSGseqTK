/*
 ============================================================================
 Name        : hello.c
 Author      : Qi Zheng
 Version     :
 Copyright   : GPL v3.0 Copyright (C) 2018  Qi Zheng
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/libsds/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */

#include "libsds.h"
#include "stdio.h"

int main(void) {
	printf("Hello World\n");
}
