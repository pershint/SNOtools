ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
COPARSER := ConfigParser.cc
ERROR := Exceptions.cc
SACPROGRAM := Sacrifice.cc
BIFURPROGRAM := Bifurcator.cc
SACOUT := Sacrifice
BIFUROUT := Bifurcator

standard:
	g++ $(ROOT_CONFIG) lib/$(SACPROGRAM)  lib/$(COPARSER) -o bin/$(SACOUT) $(ROOT_LIBS)
	g++ $(ROOT_CONFIG) lib/$(BIFURPROGRAM)  lib/$(COPARSER) -o bin/$(BIFUROUT) $(ROOT_LIBS) 
