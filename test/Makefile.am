AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/util/UCSC \
-I$(top_srcdir)/src/util/GFF \
-I$(top_srcdir)/src/util/SAMtools \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CXXFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/util/UCSC \
-I$(top_srcdir)/src/util/GFF \
-I$(top_srcdir)/src/util/SAMtools \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include \
-I$(top_srcdir)/src/libsds/include \
$(BOOST_CPPFLAGS) \
$(OPENMP_CFLAGS)

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS) $(OPENMP_CXXFLAGS)

check_PROGRAMS = \
DNAseq_test \
PrimarySeq_test \
SeqIO_test \
FMIndex_test \
Genome_test \
MetaGenomeAnno_test \
MEM_test \
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

FMIndex_test_SOURCES = FMIndex_test.cpp
FMIndex_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.la \
$(top_srcdir)/src/util/libEGUtil.a

Genome_test_SOURCES = Genome_test.cpp $(top_srcdir)/src/MSGseqTKEnv.cpp
Genome_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a \
$(BOOST_REGEX_LIB)

MetaGenomeAnno_test_SOURCES = MetaGenomeAnno_test.cpp $(top_srcdir)/src/MSGseqTKEnv.cpp
MetaGenomeAnno_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a \
$(BOOST_REGEX_LIB)

MEM_test_SOURCES = MEM_test.cpp 
MEM_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.la \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a \
$(BOOST_REGEX_LIB)

Alignment_test_SOURCES = Alignment_test.cpp 
Alignment_test_LDADD = \
$(top_srcdir)/src/libMSGseqTK.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libsds/src/libsds.la \
$(top_srcdir)/src/util/libEGUtil.a \
$(top_srcdir)/src/util/UCSC/libEGUCSC.a \
$(top_srcdir)/src/util/SAMtools/libEGSAMtools.a \
$(BOOST_REGEX_LIB) \
-lhts

TESTS = \
DNAseq_test \
PrimarySeq_test \
SeqIO_test \
FMIndex_test \
Genome_test \
MetaGenomeAnno_test \
MEM_test \
Alignment_test
