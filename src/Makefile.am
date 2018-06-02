SUBDIRS = libdivsufsort libcds util math

AM_CFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libcds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS) $(OPENMP_CXXFLAGS)

# Build convenient non-install libraries
noinst_LIBRARIES = libMSGseqClean.a

libMSGseqClean_a_SOURCES = \
DNAalphabet.cpp \
DNAseq.cpp \
PrimarySeq.cpp \
SeqIO.cpp \
RFMIndex.cpp \
MetaGenome.cpp

# Build the executables
bin_PROGRAMS = \
msgseqclean \
msgseqclean-build \
msgseqclean-inspect

msgseqclean_SOURCES = msgseqclean.cpp
msgseqclean_LDADD = \
libMSGseqClean.a \
util/libEGUtil.a \
math/libEGMath.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libcds/src/libcds.la \
$(BOOST_IOSTREAMS_LIB)

msgseqclean_build_SOURCES = msgseqclean-build.cpp MSGseqCleanEnv.cpp
msgseqclean_build_LDADD = \
libMSGseqClean.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libcds/src/libcds.la \
$(BOOST_IOSTREAMS_LIB)

msgseqclean_inspect_SOURCES = msgseqclean-inspect.cpp MSGseqCleanEnv.cpp
msgseqclean_inspect_LDADD = \
libMSGseqClean.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libcds/src/libcds.la
