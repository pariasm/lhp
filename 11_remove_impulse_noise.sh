#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

echo "	This is a SECONDARY script called by 10_preprocess_noise.h"
