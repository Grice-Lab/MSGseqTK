AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/UCSC \
-I$(top_srcdir)/SAMtools \
$(HTSLIB_CPPFLAGS)

AM_CXXFLAGS = \
-I$(top_srcdir) \
-I$(top_srcdir)/UCSC \
-I$(top_srcdir)/SAMtools \
$(HTSLIB_CPPFLAGS)

AM_CPPFLAGS = \
$(BOOST_CPPFLAGS) \
$(HTSLIB_CPPFLAGS)

check_PROGRAMS = \
GFF_test \
SAMtools_test

GFF_test_SOURCES = GFF_test.cpp
GFF_test_LDADD = \
$(top_srcdir)/UCSC/libEGUCSC.a \
$(top_srcdir)/libEGUtil.a

SAMtools_test_SOURCES = SAMtools_test.cpp
SAMtools_test_LDFLAGS = \
$(HTSLIB_LDFLAGS)

SAMtools_test_LDADD = \
$(top_srcdir)/libEGUtil.a \
$(top_srcdir)/SAMtools/libEGSAMtools.a \
-lhts

TESTS = \
GFF_test \
SAMtools_test
