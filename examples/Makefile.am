AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -I$(top_srcdir)/include \
$(OPENMP_CFLAGS)

AM_CPPFLAGS = \
$(OPENMP_CFLAGS)

noinst_PROGRAMS = bwt mksary sasearch suftest unbwt

if BUILD64
PKG_NAME=libdivsufsort64
else
PKG_NAME=libdivsufsort
endif

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

