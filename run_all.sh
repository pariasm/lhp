#! /bin/bash

# check correct number of input args
if [ "$#" -ne 1 ]; 
then
	echo "Usage: $0 sequence-folder" >&2
	exit 1
fi

SEQUENCE=$1

echo "Running full pipeline for sequence $SEQUENCE"


./10_preprocess_noise.sh $SEQUENCE
./20_stabilize_video.sh $SEQUENCE
./30_compute_optical_flow.sh $SEQUENCE
./40_run_denoising.sh $SEQUENCE
