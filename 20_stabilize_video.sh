#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Stabilizing sequence $SEQUENCE. Output stored in output_data/2_stabilization/$SEQUENCE/"


# for now, we just create sym links in the output folder
INPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
OUTPUT_DIR="output_data/2_stabilization/$SEQUENCE"
mkdir -p $OUTPUT_DIR
BASE_DIR=$(pwd)
cd $OUTPUT_DIR
for i in $(ls ${BASE_DIR}/${INPUT_DIR}/*tif);
do
	ln -s $i
done
cd -
