lib_LIBRARIES = libEGSAMtools.a

# Sources for libEGSAMtools.a
libEGSAMtools_a_SOURCES = \
BAMheader.cpp \
HTSindex.cpp \
BAM.cpp \
SAMfile.cpp

libEGSAMtools_a_CPPFLAGS = \
$(BOOST_CPPFLAGS) \
$(HTSLIB_CPPFLAGS)

libEGSAMtools_a_LIBADD = \
-lhts