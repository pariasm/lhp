#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1
echo "Preprocessing sequence $SEQUENCE. Output stored in output_data/1_preprocessing/$SEQUENCE/"

# do things e.g. call helper script
./11_remove_impulse_noise.sh $SEQUENCE
