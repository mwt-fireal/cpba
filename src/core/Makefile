# Makefile for CPBA Toolset 1.1c (CORE LIBRARY)

#we use c++
CC=g++

PLATFORM=_POSIX
UNIXPLATFORM=$(shell uname -s)

SunOS_SW_ROOT=/export/sw
Darwin_SW_ROOT=/opt
Linux_SW_ROOT=/opt

SW_ROOT=$($(UNIXPLATFORM)_SW_ROOT)

BOOSTLIB=$(SW_ROOT)/boost_1_48_0/stage/lib
LDFLAGS=-L$(SW_ROOT)/gsl-1.15/lib -L$(SW_ROOT)/boost_1_48_0/stage/lib
CPPFLAGS=-I$(SW_ROOT)/gsl-1.15/include -I$(SW_ROOT)/boost_1_48_0 
LIBS=$(SW_ROOT)/gsl-1.15/lib/libgsl.a $(SW_ROOT)/gsl-1.15/lib/libgslcblas.a $(BOOSTLIB)/libboost_filesystem.a $(BOOSTLIB)/libboost_system.a $(BOOSTLIB)/libboost_program_options.a -lstdc++

ifeq ($(strip $(UNIXPLATFORM)),Darwin)
	CFLAGS=-mmacosx-version-min=10.5 -c -D$(PLATFORM)
endif

ifeq ($(strip $(UNIXPLATFORM)),Linux)
	CFLAGS=-c -D$(PLATFORM)
endif

ifeq ($(strip $(UNIXPLATFORM)),SunOS)
	CFLAGS=-c -D$(PLATFORM)
endif
	
SOURCES=$(shell ls *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) library

.cpp.o:
	$(CC) -g $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) $< -o $@

library: $(OBJECTS)
	ar rcs libcpba.a *.o
	
# build the correct shared library
ifeq ($(strip $(UNIXPLATFORM)),Darwin)
	$(CC) -dynamiclib -w $(LIBS) $(LDFLAGS) $(OBJECTS) -o libcpba.so
endif

ifeq ($(strip $(UNIXPLATFORM)),SunOS)
	/usr/ccs/bin/ld -G -L/usr/local/lib $(LIBS) $(LDFLAGS) $(OBJECTS) -o libcpba.so
endif

clean:
	rm -rf $(OBJECTS) libcpba.a libcpba.so
