AC_PREREQ([2.69])
AC_INIT([libdivsufsort],[2.0.0],[e00011027@gmail.com])
AC_CONFIG_SRCDIR([lib/utils.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_INIT_AUTOMAKE

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_TYPE_INT64_T
AC_C_INLINE

# Checks for library functions.
#AC_FUNC_MALLOC # ignored to avoid cross-compilation problems
AC_CHECK_FUNCS([malloc])
AC_PROG_RANLIB
LT_INIT

AC_CONFIG_FILES([Makefile lib/Makefile examples/Makefile])

AC_OPENMP

AC_OUTPUT
