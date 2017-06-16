#! /bin/bash

# zoom factor
ZF=2

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

UPSA="src/utils/imscript/bin/upsa"
PLAMBDA="src/utils/imscript/bin/plambda"

SEQUENCE=$1
INPUT_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
OUTPUT_DIR="output_data/3_oflow/$SEQUENCE"
echo "	Upscaling sequence $INPUT_DIR. Output stored in $OUTPUT_DIR"

# for now, we just create sym links in the output folder
mkdir -p $OUTPUT_DIR
for i in $(ls ${INPUT_DIR}/*flo);
do
	$UPSA $ZF 2 $i $OUTPUT_DIR/$(basename $i)
	$PLAMBDA $OUTPUT_DIR/$(basename $i) "x $ZF *" -o $OUTPUT_DIR/$(basename $i)
done

