#! /bin/bash

# zoom factor
ZF=2

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

DOWNSA="src/utils/imscript/bin/downsa"


SEQUENCE=$1
OUTPUT_DIR="output_data/3_oflow/$SEQUENCE/downscaled/"
INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
echo "	Downscaling sequence $INPUT_DIR. Output stored in $OUTPUT_DIR"

# for now, we just create sym links in the output folder
mkdir -p $OUTPUT_DIR
BASE_DIR=$(pwd)
for i in $(ls ${BASE_DIR}/${INPUT_DIR}/*tif);
do
	$DOWNSA v $ZF $i $OUTPUT_DIR/$(basename $i)
done


