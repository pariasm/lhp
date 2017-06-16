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
echo "Denoising sequence $SEQUENCE. Output stored output_data/4_denoising/$SEQUENCE/"

INPUT_DIR="output_data/2_stabilization/$SEQUENCE"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/*.tif | wc -l)
	L=$((F + N - 1))
fi
