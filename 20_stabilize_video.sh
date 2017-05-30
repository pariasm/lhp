#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Stabilizing sequence $SEQUENCE. Output stored in output_data/2_stabilization/$SEQUENCE/"
