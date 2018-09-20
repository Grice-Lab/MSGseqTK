dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([libsds], [1.0], [e00011027@gmail.com])

AC_CANONICAL_SYSTEM

AC_CONFIG_MACRO_DIR([m4])

dnl Initialize automake
AM_INIT_AUTOMAKE

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(11)

dnl Initialize Libtool
LT_INIT

AC_CONFIG_FILES(Makefile
                test/Makefile
                src/Makefile
                include/Makefile
                examples/Makefile)
AC_OUTPUT