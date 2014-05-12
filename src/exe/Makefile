# Makefile for CPBA Toolset 1.1c (EXECUTBALES)

#we use c++
CC=g++

PLATFORM=_POSIX
UNIXPLATFORM=$(shell uname -s)

SunOS_SW_ROOT=/export/sw
Darwin_SW_ROOT=/opt
Linux_SW_ROOT=/opt

SW_ROOT=$($(UNIXPLATFORM)_SW_ROOT)

BOOSTLIB=$(SW_ROOT)/boost_1_48_0/stage/lib
LDFLAGS=-L$(SW_ROOT)/gsl-1.15/lib -L$(SW_ROOT)/boost_1_48_0/stage/lib -L../core
CPPFLAGS=-I$(SW_ROOT)/gsl-1.15/include -I$(SW_ROOT)/boost_1_48_0 -I../core
LIBS=$(SW_ROOT)/gsl-1.15/lib/libgsl.a $(SW_ROOT)/gsl-1.15/lib/libgslcblas.a $(BOOSTLIB)/libboost_filesystem.a $(BOOSTLIB)/libboost_system.a $(BOOSTLIB)/libboost_program_options.a -lstdc++

ifeq ($(strip $(UNIXPLATFORM)),Darwin)
	CFLAGS=-mmacosx-version-min=10.5 -D$(PLATFORM)
	EXTRALIBS=
endif

ifeq ($(strip $(UNIXPLATFORM)),Linux)
	CFLAGS=-D$(PLATFORM)
endif

ifeq ($(strip $(UNIXPLATFORM)),SunOS)
	CFLAGS=-D$(PLATFORM)
	EXTRALIBS=-L/usr/local/lib
endif
	

SOURCES=$(shell ls *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS)

.cpp.o:

	$(CC) -g $(CFLAGS) $(EXTRALIBS) $(LDFLAGS) $(CPPFLAGS) ../core/*.o $< -o $* $(LIBS)

clean:
	rm -rf $(OBJECTS) main dose_scale 
