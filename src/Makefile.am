SUBDIRS = util

AM_CFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util
#-I$(top_srcdir)/src/math \
#-I$(top_srcdir)/src/libdivsufsort/include \
#-I$(top_srcdir)/src/libcds/include \
#$(BOOST_CPPFLAGS)

AM_CXXFLAGS = -I$(top_srcdir)/src \
-I$(top_srcdir)/src/util
#-I$(top_srcdir)/src/math \
#-I$(top_srcdir)/src/libdivsufsort/include \
#-I$(top_srcdir)/src/libcds/include \
#$(BOOST_CPPFLAGS)

#AM_CPPFLAGS = $(BOOST_CPPFLAGS)
#AM_LDFLAGS = $(BOOST_LDFLAGS)

# Build the executables
bin_PROGRAMS = \
msgseqclean

msgseqclean_SOURCES = msgseqclean.cpp
#msgseqclean_CPPFLAGS = $(OPENMP_CFLAGS) $(BOOST_CPPFLAGS)
#msgseqclean_CXXFLAGS = $(OPENMP_CFLAGS)
msgseqclean_LDADD = util/libEGUtil.a
#math/libEGMath.a \
#libdivsufsort/lib/libdivsufsort.a libcds/src/libcds.la \
#$(BOOST_IOSTREAMS_LIB)

#msgseqclean_LDFLAGS = $(BOOST_LDFLAGS)