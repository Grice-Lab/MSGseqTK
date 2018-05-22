dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT(MSGseqClean, 1.0)
LT_INIT([shared static])

AC_CANONICAL_SYSTEM
AC_CANONICAL_HOST

AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT8_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_PROG_LIBTOOL
LT_INIT
#AC_PROG_RANLIB
AC_CHECK_HEADERS
AC_CHECK_HEADER_STDBOOL

AM_INIT_AUTOMAKE()

AC_LANG(C++)
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_CONFIG_MACRO_DIR([m4])

AX_CXX_COMPILE_STDCXX(11)

#AC_CHECK_FUNCS([sqrt isnan floor])
#AC_CHECK_FUNCS([strerror])

#AX_CHECK_ZLIB()
#AX_BOOST_BASE([1.50.0])
#AX_BOOST_IOSTREAMS()

AC_CONFIG_FILES([Makefile
				 src/Makefile
				 test/Makefile])

AC_CONFIG_SUBDIRS([src/util
                   src/math
                   src/libcds
                   src/libdivsufsort])

# customized macros
AC_ARG_ENABLE([build64], AS_HELP_STRING([--disable-build64], [do not build a 64 bits support]),
	[AC_DEFINE([BUILD_DIVSUFSORT64], [0], [build a 32 bits version of libdivsufsort])],
	[AC_DEFINE([BUILD_DIVSUFSORT64], [1], [build a 64 bits version of libdivsufsort])]
)

#AC_OPENMP

AC_OUTPUT
