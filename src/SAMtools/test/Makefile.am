AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/SAMtools

AM_CXXFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/SAMtools

check_PROGRAMS = \
SAMtools_test

SAMtools_test_SOURCES = SAMtools_test.cpp
SAMtools_test_LDFLAGS = \
$(HTSLIB_LDFLAGS)

SAMtools_test_LDADD = \
$(top_srcdir)/libEGSAMtools.a \
$(BOOST_REGEX_LIB) \
-lhts


TESTS = \
SAMtools_test
