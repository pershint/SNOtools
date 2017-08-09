ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
PROGRAM := BifurcatedAnalysis.cc
OUTPUT_FILE := BifurcatedAnalysis 

standard:
	g++ $(ROOT_CONFIG)  $(PROGRAM) -o $(OUTPUT_FILE) $(ROOT_LIBS)
