#-------------------------------------------------------------------------------
# This file governs the building and installation process of the TIFA library
# and its associated programs.
#
# Options given in this file follow the syntax: OPTIONS = <value>
#
# If an option is not defined in this file, its default value will be used.
# Uncomment and edit this file to suit your needs.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#                     Installation directories options
#-------------------------------------------------------------------------------

#
# Base directory serving as prefix for the other intallation directories.
#
# Default: '/usr/local'
#-------------------------------------------------------------------------------
#PREFIX = '/home/username'

#
# Install all programs in PREFIX/BINDIR if BINDIR is not an absolute path.
# If BINDIR is an absolute path, install all programs in BINDIR.
#
# Default: 'bin'
#-------------------------------------------------------------------------------
#BINDIR = 'bin'

#
# Install TIFA's C headers in PREFIX/INCLUDEDIR if INCLUDEDIR is not an
# absolute path.
# If INCLUDEDIR is an absolute path, install the headers in INCLUDEDIR.
#
# Quite a few headers are installed, so it is advised to put them all in
# a TIFA dedicated folder such as /usr/include/tifa or ~/include/tifa.
#
# Default: 'include/tifa'
#-------------------------------------------------------------------------------
#INCLUDEDIR = 'include/tifa'

#
# Install all scripts in PREFIX/SCRIPTDIR if SCRIPTDIR is not an absolute path.
# If SCRIPTDIR is an absolute path, install all scripts in SCRIPTDIR.
#
# Default: 'bin'
#-------------------------------------------------------------------------------
#SCRIPTDIR = 'bin'

#
# Install all perl module in PREFIX/PMDIR if PMDIR is not an absolute path.
# If PMDIR is an absolute path, install perl modules in PMDIR.
#
# Default: Perl's $Config{installsitelib} as defined in the Config.pm module
#-------------------------------------------------------------------------------
#PERLMODDIR = '/home/username/lib/perl5/5.8.8'

#
# Install programs or scripts' configuration files in PREFIX/CONFDIR if
# CONFDIR is not an absolute path.
# If CONFDIR is an absolute path, install all configuration files in CONFDIR.
#
# Default: 'share/tifa/conf'
#-------------------------------------------------------------------------------
#CONFDIR = 'conf'

#
# Install the TIFA library in PREFIX/LIBDIR if LIBDIR is not an absolute path.
# If LIBDIR is an absolute path, install the library in LIBDIR.
#
# Default: 'lib'
#-------------------------------------------------------------------------------
#LIBDIR = 'lib'

#
# Install TIFA's man pages in PREFIX/MANDIR if MANDIR is not an absolute path.
# If MANDIR is an absolute path, install the man pages in MANDIR.
#
# Default: 'man'
#-------------------------------------------------------------------------------
#MANDIR = 'man'

#
# Install TIFA's Doxygen documentation in PREFIX/DOCDIR if DOCDIR is not an
# absolute path.
# If DOCDIR is an absolute path, install the Doxygen documentation in DOCDIR.
#
# Default: 'share/tifa/doc'
#-------------------------------------------------------------------------------
#DOCDIR = 'doc'

#-------------------------------------------------------------------------------
#                               Build options
#-------------------------------------------------------------------------------

#
# Project wide C compiler flags.
#
# Default: '-O3'
#-------------------------------------------------------------------------------
#CCFLAGS = '-g -Wall -Wextra'

#
# Build the TIFA library and programs in the BUILDDIR directory (or in a
# subdirectory of BUILDDIR if BUILDDIR_DEPENDS_ON_PLATFORM is set to 'yes').
#
# Default: 'build'
#-------------------------------------------------------------------------------
#BUILDDIR = 'build'

#
# Set to 'yes' to build the TIFA library and programs in an platform dependent
# subdirectory of BUILDDIR.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#BUILDDIR_DEPENDS_ON_PLATFORM = 'no'

#
# Directory where the GMP library is installed. If left empty, SCons will
# try to find the GMP library in the standard locations (/usr/lib,
# /usr/local/lib, etc.)
#
# Default: '' (empty string)
#-------------------------------------------------------------------------------
#GMP_LIBDIR = '/home/username/lib'

#
# Directory where the GMP headers are installed. If left empty, SCons will
# try to find the GMP headers in the standard locations (/usr/include,
# /usr/local/include, etc.)
#
# Default: '' (empty string)
#-------------------------------------------------------------------------------
#GMP_INCDIR = '/home/username/include'

#
# Set to 'yes' to use some of GMP's "private" functions declared in the
# internal GMP headers (gmp-impl.h, longlong.h, etc.).
# Note that these headers are not installed by default, so if you did not
# compiled GMP by yourself, chances are that these headers are not available
# on your system.
#
# Default: 'no'
#-------------------------------------------------------------------------------
#USE_GMP_INTERNAL_FUNCS = 'yes'

#
# Directory where the internal GMP headers (gmp-impl.h, longlong.h, etc.) are
# installed. Obviously, this is only used if USE_GMP_INTERNAL_FUNCS is set
# to 'yes'.
#
# Default: '' (empty string)
#-------------------------------------------------------------------------------
#GMP_INTERNAL_INCDIR = '/home/username/include'

#
# Use the calloc and memset functions to allocate and/or initialize integer
# arrays to 0. Note that this works if and only if the binary representation
# of the value 0 is indeed only made up by 0 bits. This is indeed the case
# for integers on most architectures (x86, x86_64, PowerPC...).
# If in doubt, set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#USE_CALLOC_MEMSET = 'no'

#
# Native C type used to represent a string of bits.
#
# Default: 'uint64_t'
#-------------------------------------------------------------------------------
#BITSTRING_T = 'unsigned long int'

#-------------------------------------------------------------------------------
#                               Verbosity options
#-------------------------------------------------------------------------------

#
# Set to 'yes' to print critical error messages on stderr.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#PRINT_ERROR = 'no'

#
# Set to 'yes' to turn on output messages for some Tifa functions and programs.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#ALLOW_VERBOSE = 'no'

#
# Set to 'yes' to turn on output messages for the cfrac function.
# Has no effect if ALLOW_VERBOSE is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#VERBOSE_CFRAC = 'no'

#
# Set to 'yes' to turn on output messages for the fermat function.
# Has no effect if ALLOW_VERBOSE is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#VERBOSE_FERMAT = 'no'

#
# Set to 'yes' to turn on output messages for the qs function.
# Has no effect if ALLOW_VERBOSE is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#VERBOSE_QS = 'no'

#
# Set to 'yes' to turn on output messages for the siqs function.
# Has no effect if ALLOW_VERBOSE is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#VERBOSE_SIQS = 'no'

#
# Set to 'yes' to turn on output messages for the squfof function.
# Has no effect if ALLOW_VERBOSE is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#VERBOSE_SQUFOF = 'no'

#-------------------------------------------------------------------------------
#                               Timing options
#-------------------------------------------------------------------------------

#
# Set to 'yes' to turn on timing measurements for some Tifa functions and
# programs.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#ALLOW_TIMING = 'no'

#
# Set to 'yes' to turn on timing messages for the cfrac function.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#TIMING_CFRAC = 'no'

#
# Set to 'yes' to turn on timing messages for the fermat function.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#TIMING_FERMAT = 'no'

#
# Set to 'yes' to turn on timing messages for the qs function.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#TIMING_QS = 'no'

#
# Set to 'yes' to turn on timing messages for the siqs function.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#TIMING_SIQS = 'no'

#
# Set to 'yes' to turn on timing messages for the squfof function.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#TIMING_SQUFOF = 'no'

#-------------------------------------------------------------------------------
#                            Documentation options
#-------------------------------------------------------------------------------

#
# Set to yes to generate some Doxygen documentation. If set to 'no', no
# doxygen documentation will be generated, irrespective of the other
# GENERATE_DOXIGEN_* tags.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#GENERATE_DOXYGEN_DOC = 'yes'

#
# Set the paper type of the PDF Doxygen documentation. Should be one of the
# following values:
#
#    - 'a4'        (210 x 297 mm)
#    - 'a4wide'    (same as a4, but including the Latex a4wide package)
#    - 'letter'    (8.5 x 11 inches)
#    - 'legal'     (8.5 x 14 inches)
#    - 'executive' (7.25 x 10.5 inches)
#
# Default: 'a4wide'
#-------------------------------------------------------------------------------
#DOXYGEN_PAPER_TYPE = 'letter'

#
# Set to yes to generate the Doxygen documentation in HTML format.
# Has no effect if GENERATE_DOXYGEN_DOC is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#GENERATE_DOXYGEN_HTML = 'no'

#
# Set to yes to generate the Doxygen documentation in PDF format.
# Has no effect if GENERATE_DOXYGEN_DOC is set to 'no'.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#GENERATE_DOXYGEN_PDF = 'no'

#
# Set to yes to generate the Doxygen documentation in XML format.
# Has no effect if GENERATE_DOXYGEN_DOC is set to 'no'.
#
# Default: 'no'
#-------------------------------------------------------------------------------
#GENERATE_DOXYGEN_XML = 'yes'

#
# Set to yes to generate the Doxygen documentation in RTF format.
# Has no effect if GENERATE_DOXYGEN_DOC is set to 'no'.
#
# Default: 'no'
#-------------------------------------------------------------------------------
#GENERATE_DOXYGEN_RTF = 'yes'

#
# Set to yes to generate TIFA man pages.
#
# Default: 'yes'
#-------------------------------------------------------------------------------
#GENERATE_MAN_PAGES = 'no'


