#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Computing optical flow for sequence $SEQUENCE. Output stored in output_data/3_oflow/$SEQUENCE/"

# downsample to half resolution
./31_downsample_sequence.sh $SEQUENCE

# compute forwared and backward optical flow
IN_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
FLOW_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
NFRAMES=$(ls $IN_DIR/*tif | wc -l)
./32_compute_tvl1_flow.sh $IN_DIR/%03d.tif 1 $NFRAMES $FLOW_DIR


# upsample to full resolution
./33_upsample_oflow.sh $SEQUENCE
