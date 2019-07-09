###################################################################
# Copyright (c) 2019
# Author(s): Marcus D. R. Klarqvist
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
###################################################################

# Version numbers slices from the source header
LIBVER_MAJOR_SCRIPT:=`sed -n '/const int32_t DJINN_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/djinn.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const int32_t DJINN_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/djinn.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const int32_t DJINN_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/djinn.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER := $(shell echo $(LIBVER_SCRIPT))

OPTFLAGS  := -O3 -DHAVE_ZLIB -DHAVE_ZSTD -DHAVE_LZ4
CFLAGS     = -std=c99   $(OPTFLAGS) $(DEBUG_FLAGS) -g
CXXFLAGS   = -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS) -g
CXX_SOURCE = lib/frequency_model.cpp lib/pbwt.cpp lib/djinn.cpp lib/ctx_model.cpp lib/ewah_model.cpp
C_SOURCE   = 
OBJECTS    = $(C_SOURCE:.c=.o) $(CXX_SOURCE:.cpp=.o)
DEBUG_FLAGS =
DEBUG_LIBS  =
OPENSSL_PATH = 

LIBS := -lzstd -llz4

# OS X linker doesn't support -soname, and use different extension
# see : https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/100-Articles/DynamicLibraryDesignGuidelines.html
# Use LDFLAGS for passing additional flags. Avoid using LD_LIB_FLAGS.
ifneq ($(shell uname), Darwin)
SHARED_EXT   = so
LD_LIB_FLAGS = -shared '-Wl,-rpath-link,$$ORIGIN/,-rpath-link,$(PWD),-soname,libdjinn.$(SHARED_EXT)' $(LDFLAGS)
else
SHARED_EXT   = dylib
LD_LIB_FLAGS = -dynamiclib -install_name "@rpath/libdjinn.$(SHARED_EXT)" '-Wl,-rpath,@loader_path/,-rpath,$(PWD)' $(LDFLAGS) 
endif


# Default target
all: djinn

debug_size: DEBUG_FLAGS += -DDEBUG_SIZE
debug_size: djinn
debug: DEBUG_FLAGS += -DDEBUG_PBWT -DDEBUG_WAH -DDEBUG_CONTEXT -DDEBUG_SIZE -g
debug: DEBUG_LIBS += -lcrypto
debug: djinn

# Generic rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -fPIC -c -DVERSION=\"$(GIT_VERSION)\" -o $@ $<

library: $(OBJECTS)
	@echo 'Building dynamic library...'
	$(CXX) $(LD_LIB_FLAGS) $(LIBRARY_PATHS) $(OBJECTS) -pthread $(LIBS) -o libdjinn.$(SHARED_EXT).$(LIBVER)
	@echo 'Building static library...'
	$(AR) crs libdjinn.a $(OBJECTS)
	@echo 'Symlinking library...'
	ln -sf libdjinn.$(SHARED_EXT).$(LIBVER) libdjinn.$(SHARED_EXT)
	ln -sf libdjinn.$(SHARED_EXT).$(LIBVER) libdjinn.$(SHARED_EXT)

djinn: library
	$(CXX) $(CXXFLAGS) main.cpp -Ilib/ -L$(PWD) -pthread $(LIBS) -lhts -ldjinn '-Wl,-rpath,$$ORIGIN/,-rpath,$(PWD)' -o djinn

clean:
	rm -f *.o lib/*.o *.a *.$(SHARED_EXT).* *.$(SHARED_EXT) djinn

.PHONY: all library clean debug debug_size
