#! /bin/bash

# zoom factor
ZF=2

# check correct number of input args
if [ "$#" -lt 1 ];
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

F=${2:-1} # first frame
L=${3:-0} # last frame

SEQUENCE=$1
OUTPUT_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
echo "	Downscaling sequence $INPUT_DIR. Output stored in $OUTPUT_DIR"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# downsample the sequence
DOWNSA="src/utils/imscript/bin/downsa"
mkdir -p $OUTPUT_DIR
BASE_DIR=$(pwd)
for i in $(seq $F $L)
do
	N=$(printf %03d $i)
	$DOWNSA v $ZF ${INPUT_DIR}/$N.tif $OUTPUT_DIR/$N.tif
done


