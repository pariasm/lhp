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
INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
OFLOW_DIR="output_data/3_oflow_dw0.10/$SEQUENCE"
OUTPUT_DIR="output_data/4_denoising/$SEQUENCE/rbilf-dw0.10"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

echo "Denoising sequence $SEQUENCE, from frame $F to $L."
echo "Output stored $OUTPUT_DIR"

# create directory
mkdir -p $OUTPUT_DIR

# nldct binary
DENO="src/4_denoising/rbilf/build/bin/rbilf"
SIGMA=$(cat "$INPUT_DIR/sigma.txt")

# run denoising (first step)
$DENO \
-i ${INPUT_DIR}/%03d.tif -f $F -l $L -s $SIGMA -o ${OFLOW_DIR}/%03d.b.flo \
-d ${OUTPUT_DIR}/d_%03d.tif -v

