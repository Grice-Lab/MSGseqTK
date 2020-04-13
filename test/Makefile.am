AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/SAMtools \
$(HTSLIB_CPPFLAGS)

AM_CXXFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/SAMtools \
$(HTSLIB_CPPFLAGS)

AM_CPPFLAGS = \
$(BOOST_CPPFLAGS) \
$(HTSLIB_CPPFLAGS)


check_PROGRAMS = \
SAMtools_test

SAMtools_test_SOURCES = SAMtools_test.cpp

SAMtools_test_LDFLAGS = \
$(HTSLIB_LDFLAGS)

SAMtools_test_LDADD = \
$(top_srcdir)/libEGSAMtools.a \
-lhts


TESTS = \
SAMtools_test
