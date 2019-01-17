AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/include \
$(OPENMP_CFLAGS)

AM_CFLAGS = -I$(top_srcdir)/include \
$(OPENMP_CFLAGS)

if BUILD64
lib_LIBRARIES = libdivsufsort64.a
libdivsufsort64_a_SOURCES = divsufsort.c sssort.c trsort.c utils.c
else
lib_LIBRARIES = libdivsufsort.a
libdivsufsort_a_SOURCES = divsufsort.c sssort.c trsort.c utils.c
endif