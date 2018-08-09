# This is the Makefile for the EGUtil library
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS = -I$(top_srcdir)

AM_CXXFLAGS = -I$(top_srcdir)

# Build a static/shared library
lib_LIBRARIES = libEGUCSC.a
libEGUCSC_a_SOURCES = \
GFF.cpp

libEGUCSC_a_LIBADD = \
$(top_srcdir)/libEGUtil.a
