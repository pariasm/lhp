#! /bin/bash

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

F=${2:-1} # first frame
L=${3:-0} # last frame 

SEQUENCE=$1
echo "Computing optical flow for sequence $SEQUENCE. Output stored in output_data/3_oflow/$SEQUENCE/"

INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
FLOW_DIR="output_data/3_oflow/$SEQUENCE"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# compute forwared and backward optical flow
./32_compute_tvl1_flow.sh $INPUT_DIR/%03d.tif $F $L $FLOW_DIR

## the downsampling and upsampling is now done in the tvl1code,
## thus we don't need this:
## # downsample to half resolution
## ./31_downsample_sequence.sh $SEQUENCE $F $L
## 
## # compute forwared and backward optical flow
## IN_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
## FLOW_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
## ./32_compute_tvl1_flow.sh $IN_DIR/%03d.tif $F $L $FLOW_DIR
## 
## 
## # upsample to full resolution
## ./33_upsample_oflow.sh $SEQUENCE $F $L
## 
## 
## # clean intermediate data
## rm -rf $FLOW_DIR
