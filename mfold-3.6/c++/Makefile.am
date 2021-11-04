AM_CPPFLAGS = -DPKGDATADIR=\"$(pkgdatadir)/\" 
AM_CXXFLAGS = -Wall -O3 -Wno-deprecated -fomit-frame-pointer 

EXTRA_DIST = algorithm.cpp algorithm.h defines.h interface.h \
structure.cpp structure.h 

bin_PROGRAMS = efn2


efn2_SOURCES = efn2.cpp 
efn2_CXXFLAGS = $(AM_CXXFLAGS)
