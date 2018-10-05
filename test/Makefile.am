AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = 

AM_CXXFLAGS = \
-I$(top_srcdir)/include 

check_PROGRAMS = BitStr_test

BitStr_test_SOURCES = BitStr_test.cpp
BitStr_test_LDADD = \
$(top_srcdir)/src/libsds.la 

TESTS = BitStr_test