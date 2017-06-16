#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Stabilizing sequence $SEQUENCE. Output stored in output_data/2_stabilization/$SEQUENCE/"

STABI="src/2_stabilization/estadeo_1.1/bin/estadeo"

INPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
OUTPUT_DIR="output_data/2_stabilization/$SEQUENCE"
mkdir -p $OUTPUT_DIR

$STABI $INPUT_DIR/%03d.tif 1 300 -1 -1 -1 -o $OUTPUT_DIR/%03d.tif 
