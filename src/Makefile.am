SUBDIRS = libdivsufsort libsds util SAMtools math

AM_CFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/util/UCSC \
-I$(top_srcdir)/src/util/GFF \
-I$(top_srcdir)/src/SAMtools \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CXXFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/util/UCSC \
-I$(top_srcdir)/src/util/GFF \
-I$(top_srcdir)/src/SAMtools \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CXXFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CFLAGS) $(OPENMP_CXXFLAGS) $(HTSLIB_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

# Build convenient non-install libraries
noinst_LIBRARIES = libMSGseqTK.a

libMSGseqTK_a_SOURCES = \
MSGseqTK_main.cpp \
DNAalphabet.cpp \
DNAseq.cpp \
QualStr.cpp \
PrimarySeq.cpp \
SeqIO.cpp \
Loc.cpp \
GLoc.cpp \
SMEM.cpp \
SeedPair.cpp \
SeedChain.cpp \
ScoreScheme.cpp \
PairingScheme.cpp \
FMDIndex.cpp \
Genome.cpp \
MetaGenome.cpp \
MetaGenomeAnno.cpp \
Alignment.cpp

# Build the executables
bin_PROGRAMS = \
msgseqtk-build \
msgseqtk-anno \
msgseqtk-inspect \
msgseqtk-mergedb \
msgseqtk-clean \
msgseqtk-align

msgseqtk_build_SOURCES = msgseqtk-build.cpp MSGseqTKEnv.cpp
msgseqtk_build_LDADD = \
libMSGseqTK.a \
util/UCSC/libEGUCSC.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.a \
$(BOOST_IOSTREAMS_LIB) \
$(OPENMP_CXXFLAGS)

msgseqtk_anno_SOURCES = msgseqtk-anno.cpp MSGseqTKEnv.cpp
msgseqtk_anno_LDADD = \
libMSGseqTK.a \
util/UCSC/libEGUCSC.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.a \
$(BOOST_IOSTREAMS_LIB)

msgseqtk_inspect_SOURCES = msgseqtk-inspect.cpp MSGseqTKEnv.cpp
msgseqtk_inspect_LDADD = \
libMSGseqTK.a \
util/UCSC/libEGUCSC.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.a \
$(BOOST_IOSTREAMS_LIB)

msgseqtk_mergedb_SOURCES = msgseqtk-mergedb.cpp MSGseqTKEnv.cpp
msgseqtk_mergedb_LDADD = \
libMSGseqTK.a \
util/UCSC/libEGUCSC.a \
util/libEGUtil.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.a \
$(OPENMP_CXXFLAGS)

msgseqtk_clean_SOURCES = msgseqtk-clean.cpp MSGseqTKEnv.cpp
msgseqtk_clean_LDADD = \
libMSGseqTK.a \
util/UCSC/libEGUCSC.a \
util/libEGUtil.a \
math/libEGMath.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.a \
$(BOOST_IOSTREAMS_LIB) \
$(OPENMP_CXXFLAGS)

msgseqtk_align_SOURCES = msgseqtk-align.cpp MSGseqTKEnv.cpp
msgseqtk_align_LDADD = \
libMSGseqTK.a \
util/UCSC/libEGUCSC.a \
util/libEGUtil.a \
SAMtools/libEGSAMtools.a \
math/libEGMath.a \
libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
libsds/src/libsds.a \
$(BOOST_IOSTREAMS_LIB) \
$(OPENMP_CXXFLAGS) \
-lhts

msgseqtk_align_LDFLAGS = \
$(AM_LDFLAGS) \
$(HTSLIB_LDFLAGS)
