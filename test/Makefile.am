AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/util/UCSC \
-I$(top_srcdir)/src/SAMtools \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CXXFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/util/UCSC \
-I$(top_srcdir)/src/SAMtools \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CXXFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CFLAGS) $(OPENMP_CXXFLAGS) $(HTSLIB_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

check_PROGRAMS = \
DNAseq_test \
PrimarySeq_test \
SeqIO_test \
FMDIndex_test \
Genome_test \
SMEM_test \
SeedChain_test \
Alignment_test

LIBDIVSUFSORT_PKG_NAME=libdivsufsort64

DNAseq_test_SOURCES = DNAseq_test.cpp
DNAseq_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a

PrimarySeq_test_SOURCES = PrimarySeq_test.cpp
PrimarySeq_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a

SeqIO_test_SOURCES = SeqIO_test.cpp
SeqIO_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/util/libEGUtil.a

FMDIndex_test_SOURCES = FMDIndex_test.cpp
FMDIndex_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.a \
$(top_srcdir)/src/util/libEGUtil.a

Genome_test_SOURCES = Genome_test.cpp $(top_srcdir)/src/MSGseqTKEnv.cpp
Genome_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a

SMEM_test_SOURCES = SMEM_test.cpp $(top_srcdir)/src/MSGseqTKEnv.cpp
SMEM_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a

SeedChain_test_SOURCES = SeedChain_test.cpp $(top_srcdir)/src/MSGseqTKEnv.cpp
SeedChain_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a

Alignment_test_SOURCES = Alignment_test.cpp 
Alignment_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a \
$(top_srcdir)/src/SAMtools/libEGSAMtools.a \
-lhts

Alignment_test_LDFLAGS = \
$(AM_LDFLAGS) \
$(HTSLIB_LDFLAGS)

TESTS = \
DNAseq_test \
PrimarySeq_test \
SeqIO_test \
FMDIndex_test \
Genome_test \
SMEM_test \
SeedChain_test \
Alignment_test
