dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([libsds], [1.0], [e00011027@gmail.com])

AC_CANONICAL_SYSTEM

AC_CONFIG_MACRO_DIR([m4])

AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_CHECK_FUNCS([strerror])

dnl Initialize automake
AM_INIT_AUTOMAKE

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(11)
AC_PROG_RANLIB

dnl Initialize Libtool
LT_INIT

AC_CONFIG_FILES(Makefile
                test/Makefile
                src/Makefile
                include/Makefile
                examples/Makefile)
AC_OUTPUT
