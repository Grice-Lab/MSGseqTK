SUBDIRS = libdivsufsort libsds util math

AM_CFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS) \
$(OPENMP_CXXFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS) $(OPENMP_CXXFLAGS) $(BOOST_REGEX_LIB)

# Build convenient non-install libraries
noinst_LIBRARIES = libMSGseqTK.a

libMSGseqTK_a_SOURCES = \
DNAalphabet.cpp \
DNAseq.cpp \
QualStr.cpp \
PrimarySeq.cpp \
SeqIO.cpp \
Loc.cpp \
MEM.cpp \
MEMS.cpp \
FMIndex.cpp \
Genome.cpp \
MetaGenome.cpp \
GenomeAnno.cpp \
MetaGenomeAnno.cpp

# Build the executables
bin_PROGRAMS = \
msgseqtk-clean \
msgseqtk-build \
msgseqtk-inspect \
msgseqtk-mergedb

msgseqtk_clean_SOURCES = msgseqtk-clean.cpp MSGseqTKEnv.cpp
msgseqtk_clean_LDADD = \
libMSGseqTK.a \
util/libEGUtil.a \
math/libEGMath.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.la \
$(BOOST_IOSTREAMS_LIB)

msgseqtk_build_SOURCES = msgseqtk-build.cpp MSGseqTKEnv.cpp
msgseqtk_build_LDADD = \
libMSGseqTK.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.la \
$(BOOST_IOSTREAMS_LIB)

msgseqtk_inspect_SOURCES = msgseqtk-inspect.cpp MSGseqTKEnv.cpp
msgseqtk_inspect_LDADD = \
libMSGseqTK.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.la

msgseqtk_mergedb_SOURCES = msgseqtk-mergedb.cpp MSGseqTKEnv.cpp
msgseqtk_mergedb_LDADD = \
libMSGseqTK.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.la
