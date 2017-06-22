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
echo "Denoising sequence $SEQUENCE. Output stored output_data/4_denoising/$SEQUENCE/nldct/"

INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
OFLOW_DIR="output_data/3_oflow/$SEQUENCE"
OUTPUT_DIR="output_data/4_denoising/$SEQUENCE/nldct"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# nldct binary
DENO="src/4_denoising/nldct/build/bin/vnlbayes"

# create directory
mkdir -p $OUTPUT_DIR

# run denoising (first step)
$DENO \
-i ${INPUT_DIR}/%03d.tif -f $F -l $L -sigma 10 -has-noise \
-px2 0 -px1 3 -pt1 2 -wx1 54 -wt1 6 -np1 120 \
-bsic ${OUTPUT_DIR}/b_%03d.tif

# run denoising (second step)
$DENO \
-i ${INPUT_DIR}/%03d.tif -f $F -l $L -sigma 1 -has-noise \
-b ${OUTPUT_DIR}/b_%03d.tif -f $F -l $L -sigma 1 -has-noise \
-px1 0 -px2 3 -pt2 2 -wx2 54 -wt2 6 -np2 60 \
-deno ${OUTPUT_DIR}/d_%03d.tif

# clean
rm {bsic,diff}_???.png
mv measures.txt ${OUTPUT_DIR}/measures


