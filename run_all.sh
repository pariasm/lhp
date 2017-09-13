#! /bin/bash

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

F=${2:-1} # first frame
L=${3:-0} # last frame 

SEQUENCE=$1

echo "Running full pipeline for sequence $SEQUENCE"

export PATH=`pwd`/bin/:$PATH


./10_preprocess_noise.sh $SEQUENCE $F $L
./20_stabilize_video.sh $SEQUENCE $F $L
./30_compute_optical_flow.sh $SEQUENCE $F $L
./40_run_denoising.sh $SEQUENCE $F $L
./50_run_deblurring.sh $SEQUENCE $F $L

