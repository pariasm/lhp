# makefile for running the "10_preprocessing" pipeline
# this file is copied into the directory containing the links to the input files
# and then it is run by "make -f preprocess.mk"

INPUTS   = $(shell ls i*tif)
LINEAR   = $(addprefix lin_,$(INPUTS))

all: avglin.tif stdlin.tif

lin_i%.tif: i%.tif
	plambda $^ "1 + log" -o $@

avglin.tif: $(LINEAR)
	veco avg $^ -o $@

stdlin.tif: $(LINEAR)
	veco std $^ -o $@
