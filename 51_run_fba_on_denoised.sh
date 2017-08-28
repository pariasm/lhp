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
echo "Denoising sequence $SEQUENCE. Output stored output_data/4_denoising/$SEQUENCE/nlbayes/"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# fba binary
FBA="src/5_deblurring/fba/fba"

# run deblurring on denoised images
INPUT_DIR="output_data/4_denoising/$SEQUENCE/nlbayes_cbsic_l1"
OUTPUT_DIR="output_data/5_deblurring/$SEQUENCE/"
mkdir -p $OUTPUT_DIR
$FBA 3 128 11 3 1 1 ${INPUT_DIR}/d_*.tif /dev/null ${OUTPUT_DIR}/%03d.tif /dev/null

