#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Preprocessing sequence $SEQUENCE. Output stored in output_data/1_preprocessing/$SEQUENCE/"


# for now, we just create sym links in the output folder
OUTPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
mkdir -p $OUTPUT_DIR
BASE_DIR=$(pwd)
cd $OUTPUT_DIR
for i in $(ls ${BASE_DIR}/input_data/${SEQUENCE}/*tif);
do
	ln -s $i
done
cd -
