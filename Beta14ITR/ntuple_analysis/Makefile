ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
PROGRAM := B14ITR.cc
OUTPUT_FILE := B14ITR 

standard:
	g++ $(ROOT_CONFIG)  $(PROGRAM) -o $(OUTPUT_FILE) $(ROOT_LIBS)
