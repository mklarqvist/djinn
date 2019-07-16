#AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
AM_CXXFLAGS = -fPIC
if ZSTD_AVAIL
AM_CXXFLAGS += -DHAVE_ZSTD
endif
if LZ4_AVAIL
AM_CXXFLAGS += -DHAVE_LZ4
endif
if ZLIB_AVAIL
AM_CXXFLAGS += -DHAVE_ZLIB
endif
if HTS_AVAIL
AM_CXXFLAGS += -DHAVE_HTS
endif
#SUBDIRS = lib

AUTOMAKE_OPTIONS = subdir-objects
#ACLOCAL_AMFLAGS = $(ACLOCAL_FLAGS)

bin_PROGRAMS = djinn

djinn_SOURCES = main.cpp
djinn_LDADD = libdjinn.la
djinn_CXXFLAGS = -Ilib/

lib_LTLIBRARIES = libdjinn.la
libdjinn_la_LDFLAGS = -version-info 0:1:0
libdjinn_la_SOURCES = lib/compressors.h lib/ctx_model.cpp lib/djinn.cpp lib/djinn.h lib/ewah_model.cpp lib/frequency_model.cpp lib/frequency_model.h lib/pbwt.cpp lib/pbwt.h lib/vcf_reader.h