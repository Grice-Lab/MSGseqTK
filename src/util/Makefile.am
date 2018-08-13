# This is the Makefile for the EGUtil library
AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = $(BOOST_CPPFLAGS)

# Build a static/shared library
lib_LIBRARIES = libEGUtil.a
libEGUtil_a_SOURCES = \
CommandOptions.cpp \
ProgLog.cpp \
StringUtils.cpp \
TSVRecord.cpp \
TSVScanner.cpp \
VersionSequence.cpp \
ProgEnv.cpp \
UCSC/GFF.cpp

check_PROGRAMS = \
GFF_test 

GFF_test_SOURCES = $(top_srcdir)/UCSC/GFF_test.cpp
GFF_test_LDADD = $(top_srcdir)/libEGUtil.a

TESTS = \
GFF_test