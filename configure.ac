#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([DJINN], [0.1.0], [mk819@cam.ac.uk])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
: ${CXXFLAGS=-O3}


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S

# Checks for libraries.
LT_INIT
AC_ENABLE_SHARED
AC_ENABLE_STATIC
#AC_PROG_LIBTOOL(libtool)

# FIXME: Replace `main' with a function in `-ldjinn':
#AC_CHECK_LIB([djinn], [main])

zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])

# FIXME: Replace `main' with a function in `-lhts':
hts_devel=ok
AC_CHECK_HEADER([htslib/vcf.h], [], [hts_devel=missing], [;])
AC_CHECK_LIB(hts, bcf_hdr_read, [], [hts_devel=missing])

# FIXME: Replace `main' with a function in `-llz4':
lz4_devel=ok
AC_CHECK_HEADER([lz4.h], [], [lz4_devel=missing], [;])
AC_CHECK_LIB(lz4, LZ4_compress_HC,  [], [lz4_devel=missing])

# FIXME: Replace `main' with a function in `-lzstd':
zstd_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zstd.h], [], [zstd_devel=missing], [;])
AC_CHECK_LIB(zstd, ZSTD_compress,  [], [zstd_devel=missing])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([memset])

AM_CONDITIONAL([ZSTD_AVAIL], [test $zstd_devel == ok])
AM_CONDITIONAL([LZ4_AVAIL], [test $lz4_devel == ok])
AM_CONDITIONAL([HTS_AVAIL], [test $hts_devel == ok])
AM_CONDITIONAL([ZLIB_AVAIL], [test $zlib_devel == ok])

AC_MSG_NOTICE([version: $ZSTD_AVAIL.$LZ4_AVAIL.$HTS_AVAIL])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
