#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Computing optical flow for sequence $SEQUENCE. Output stored in output_data/1_oflow/$SEQUENCE/"
