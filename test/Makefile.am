AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = \
-I$(top_srcdir)/src \
-I$(top_srcdir)/src/util \
-I$(top_srcdir)/src/math \
-I$(top_srcdir)/src/libdivsufsort/include
#-I$(top_srcdir)/src/libcds/include \
#$(BOOST_CPPFLAGS) \
#$(OPENMP_CFLAGS)

#AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(OPENMP_CFLAGS)
#AM_LDFLAGS = $(BOOST_LDFLAGS) $(OPENMP_CFLAGS)

#$(top_srcdir)/src/libMSGseqClean_common.a
#$(top_srcdir)/util/libEGUtil.a \
#$(top_srcdir)/math/libEGMath.a \
#$(top_srcdir)/src/libdivsufsort/lib/libdivsufsort.a \
#$(top_srcdir)/src/libcds/src/libcds.la \
#$(BOOST_IOSTREAMS_LIB)

check_PROGRAMS = \
DNAseq_test \
PrimarySeq_test \
SeqIO_test \
RFMIndex_test

LIBDIVSUFSORT_PKG_NAME=libdivsufsort64

DNAseq_test_SOURCES = DNAseq_test.cpp
DNAseq_test_LDADD = \
$(top_srcdir)/src/libMSGseqClean_common.a

PrimarySeq_test_SOURCES = PrimarySeq_test.cpp
PrimarySeq_test_LDADD = \
$(top_srcdir)/src/libMSGseqClean_common.a

SeqIO_test_SOURCES = SeqIO_test.cpp
SeqIO_test_LDADD = \
$(top_srcdir)/src/libMSGseqClean_common.a

RFMIndex_test_SOURCES = RFMIndex_test.cpp
RFMIndex_test_LDADD = \
$(top_srcdir)/src/libMSGseqClean_common.a \
$(top_srcdir)/src/libdivsufsort/lib/${LIBDIVSUFSORT_PKG_NAME}.a \
$(top_srcdir)/src/libcds/src/libcds.la

TESTS = \
DNAseq_test \
PrimarySeq_test \
SeqIO_test \
RFMIndex_test