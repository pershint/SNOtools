ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
COPARSER := ConfigParser.cc
PROGRAM := Sacrifice.cc
TESTOUT := Sacrifice

standard:
	g++ $(ROOT_CONFIG) lib/$(PROGRAM)  lib/$(COPARSER) -o bin/$(TESTOUT) $(ROOT_LIBS) 
