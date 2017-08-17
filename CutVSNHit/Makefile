ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
PROGRAM := CutVSNHit.cc
OUTPUT_FILE := CutVSNHit

standard:
	g++ $(ROOT_CONFIG)  $(PROGRAM) -o $(OUTPUT_FILE) $(ROOT_LIBS)
