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

OPTFLAGS  := -O3 -DLZ4_AVAIL -DZSTD_AVAIL
CFLAGS     = -std=c99 $(OPTFLAGS) $(DEBUG_FLAGS) -g
CPPFLAGS   = -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS) -g
CPP_SOURCE = frequency_model.cpp main.cpp pbwt.cpp gt_compressor.cpp djinn.cpp ctx_model.cpp ewah_model.cpp
C_SOURCE   = 
OBJECTS    = $(C_SOURCE:.c=.o) $(CPP_SOURCE:.cpp=.o)
DEBUG_FLAGS =
DEBUG_LIBS  =
OPENSSL_PATH = 

# Default target
all: djinn

debug_size: DEBUG_FLAGS += -DDEBUG_SIZE
debug_size: djinn
debug: DEBUG_FLAGS += -DDEBUG_PBWT -DDEBUG_WAH -DDEBUG_CONTEXT -DDEBUG_SIZE -g
debug: DEBUG_LIBS += -lcrypto
debug: djinn

# Generic rules
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(OPENSSL_PATH) -c -o $@ $<

djinn: $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(OBJECTS) $(OPENSSL_PATH) -lzstd -llz4 -lhts $(DEBUG_LIBS) -o djinn

clean:
	rm -f $(OBJECTS)
	rm -f djinn

.PHONY: all clean debug debug_size
