AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = 

AM_CXXFLAGS = \
-I$(top_srcdir)/include 

check_PROGRAMS = \
BitStr_test \
BitSeqGGMN_test \
BitSeqRRR_test \
WaveletTreeRRR_test

BitStr_test_SOURCES = BitStr_test.cpp
BitStr_test_LDADD = \
$(top_srcdir)/src/libsds.la 

BitSeqGGMN_test_SOURCES = BitSeqGGMN_test.cpp
BitSeqGGMN_test_LDADD = \
$(top_srcdir)/src/libsds.la 

BitSeqRRR_test_SOURCES = BitSeqRRR_test.cpp
BitSeqRRR_test_LDADD = \
$(top_srcdir)/src/libsds.la 

WaveletTreeRRR_test_SOURCES = WaveletTreeRRR_test.cpp
WaveletTreeRRR_test_LDADD = \
$(top_srcdir)/src/libsds.la 

TESTS = \
BitStr_test \
BitSeqGGMN_test \
BitSeqRRR_test \
WaveletTreeRRR_test
