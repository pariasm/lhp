#! /bin/bash

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

F=${2:-1} # first frame
L=${3:-0} # last frame 

# temporal size of patch
PT=2

SEQUENCE=$1
INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
OFLOW_DIR="output_data/3_oflow/$SEQUENCE"
if [ $PT == 1 ];
then
	OUTPUT_DIR="output_data/4_denoising/$SEQUENCE/nldct-8x1"
else
	OUTPUT_DIR="output_data/4_denoising/$SEQUENCE/nldct-6x2"
fi

echo "Denoising sequence $SEQUENCE. Output stored in $OUTPUT_DIR"


# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# create directory
mkdir -p $OUTPUT_DIR

# nldct binary
DENO="src/4_denoising/nldct/build/bin/nldct"
SIGMA=$(cat "$INPUT_DIR/sigma.txt")

# run denoising (first step)
$DENO \
-i ${INPUT_DIR}/%03d.tif -f $F -l $L -sigma $SIGMA -has-noise \
-fof ${OFLOW_DIR}/%03d.f.flo -bof ${OFLOW_DIR}/%03d.b.flo \
-px2 0 -pt1 $PT -wt1 4 -bsic ${OUTPUT_DIR}/b_%03d.tif

mv measures.txt ${OUTPUT_DIR}/measures_basic

# run denoising (second step)
$DENO \
-i ${INPUT_DIR}/%03d.tif -f $F -l $L -sigma $SIGMA -has-noise \
-b ${OUTPUT_DIR}/b_%03d.tif \
-fof ${OFLOW_DIR}/%03d.f.flo -bof ${OFLOW_DIR}/%03d.b.flo \
-px1 0 -pt2 $PT -wt2 4 -deno ${OUTPUT_DIR}/d_%03d.tif

# clean
rm {bsic,diff}_???.png
mv measures.txt ${OUTPUT_DIR}/measures

