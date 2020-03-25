SUBDIRS=. test
  
# This is the Makefile for the EGUtil library
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

lib_LIBRARIES = libEGSAMtools.a

# Sources for libEGSAMtools.a
libEGSAMtools_a_SOURCES = \
BAMheader.cpp \
HTSindex.cpp \
BAM.cpp \
SAMfile.cpp
