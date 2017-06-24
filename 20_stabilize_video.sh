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
echo "Stabilizing sequence $SEQUENCE. Output stored in output_data/2_stabilization/$SEQUENCE/"

STABI="src/2_stabilization/estadeo_1.1/bin/estadeo"

INPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
OUTPUT_DIR="output_data/2_stabilization/$SEQUENCE"
mkdir -p $OUTPUT_DIR

# determine extension
FF=$(printf %03d $F)
if [ -f "input_data/${SEQUENCE}/${FF}.tif" ]
then
	EXT="tif"
elif [ -f "input_data/${SEQUENCE}/${FF}.tiff" ]
then
	EXT="tiff"
elif [ -f "input_data/${SEQUENCE}/${FF}.png" ]
then
	EXT="png"
else
	echo "File ${SEQUENCE}/${FF}.{tif,tiff,png} not found"
	exit 1
fi

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.$EXT | wc -l)
	L=$((F + N - 1))
fi

$STABI $INPUT_DIR/%03d.$EXT $F $L -1 -1 -1 -o $OUTPUT_DIR/%03d.tif 

# clean
rm $INPUT_DIR/*.tif
