ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
PROGRAM := DCSacrifice.cc
OUTPUT_FILE := DCSacrifice

standard:
	g++ $(ROOT_CONFIG)  $(PROGRAM) -o $(OUTPUT_FILE) $(ROOT_LIBS)
