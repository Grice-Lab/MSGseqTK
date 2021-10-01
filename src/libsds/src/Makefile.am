# This is the Makefile for the EGriceLab::Math library
ACLOCAL_AMFLAGS = -I m4

lib_LIBRARIES = libsds.a

#######################################
# Build information for each library

# Sources for libsds
libsds_a_SOURCES = \
BitSeq.cpp \
BitSeqGGMN.cpp \
BitSeqRRR.cpp \
Seq.cpp \
WaveletTreeRRR.cpp

# Linker options libTestProgram
libsds_a_LDFLAGS =

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libsds_a_CPPFLAGS = -I$(top_srcdir)/include
