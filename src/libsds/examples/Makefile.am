#######################################
# The list of executables we are building seperated by spaces
# the 'bin_' indicates that these build products will be installed
# in the $(bindir) directory. For example /usr/bin
#bin_PROGRAMS=exampleProgram

# Because a.out is only a sample program we don't want it to be installed.
# The 'noinst_' prefix indicates that the following targets are not to be
# installed.
noinst_PROGRAMS = \
hello \
hello2

#######################################
# Build information for each executable. The variable name is derived
# by use the name of the executable with each non alpha-numeric character is
# replaced by '_'. So a.out becomes a_out and the appropriate suffex added.
# '_SOURCES' for example.

ACLOCAL_AMFLAGS=-I ../m4

# Sources for the a.out 
hello_SOURCES= hello.c

# Libraries for a.out
hello_LDADD = $(top_srcdir)/src/libsds.a

# Linker options for a.out
hello_LDFLAGS = -rpath `cd $(top_srcdir);pwd`/src/.libs

# Compiler options for a.out
hello_CPPFLAGS = -I$(top_srcdir)/include


hello2_SOURCES = hello.cpp

hello2_LDADD = $(top_srcdir)/src/libsds.a

hello2_CPPFLAGS = -I$(top_srcdir)/include 
