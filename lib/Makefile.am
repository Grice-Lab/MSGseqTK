AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/include \
$(OPENMP_CFLAGS)

AM_CFLAGS = -I$(top_srcdir)/include \
$(OPENMP_CFLAGS)

lib_LIBRARIES = libdivsufsort.a
libdivsufsort_a_SOURCES = divsufsort.c sssort.c trsort.c utils.c
