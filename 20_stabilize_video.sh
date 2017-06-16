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

# determine extension
if [ -f "${INPUT_DIR}/001.tif" ]
then
	EXT="tif"
elif [ -f "${INPUT_DIR}/001.tiff" ]
then
	EXT="tiff"
elif [ -f "${INPUT_DIR}/001.png" ]
then
	EXT="png"
else
	echo "File ${INPUT_DIR}/001.{tif,tiff,png} not found"
	exit 1
fi

NFRAMES=$(ls $INPUT_DIR/*$EXT | wc -l)

$STABI $INPUT_DIR/%03d.$EXT 1 $NFRAMES -1 -1 -1 -o $OUTPUT_DIR/%03d.tif 
