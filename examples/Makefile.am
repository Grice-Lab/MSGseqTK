AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -I$(top_srcdir)/include \
$(OPENMP_CFLAGS)

AM_CPPFLAGS = \
$(OPENMP_CFLAGS)

check_PROGRAMS = \
bwt \
mksary \
sasearch \
suftest \
unbwt \
bwt64 \
mksary64 \
sasearch64 \
suftest64 \
unbwt64

PKG_NAME=libdivsufsort

bwt_SOURCES = bwt.c
bwt_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

mksary_SOURCES = mksary.c
mksary_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

sasearch_SOURCES = sasearch.c
sasearch_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

suftest_SOURCES = suftest.c
suftest_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

unbwt_SOURCES = unbwt.c
unbwt_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

bwt64_SOURCES = bwt64.c
bwt64_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

mksary64_SOURCES = mksary64.c
mksary64_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

sasearch64_SOURCES = sasearch64.c
sasearch64_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

suftest64_SOURCES = suftest64.c
suftest64_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

unbwt64_SOURCES = unbwt64.c
unbwt64_LDADD = $(top_srcdir)/lib/${PKG_NAME}.a

TESTS = \
bwt \
mksary \
sasearch \
suftest \
unbwt \
bwt64 \
mksary64 \
sasearch64 \
suftest64 \
unbwt64
