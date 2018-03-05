#! /bin/bash

SEQUENCE=$1
F=${2:-1} # first frame
L=${3:-0} # last frame 

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

echo "Denoising sequence $SEQUENCE. Output stored output_data/4_denoising/$SEQUENCE/"

INPUT_DIR="output_data/2_stabilization/$SEQUENCE"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

## # run ponomarenko's noise estimator
#./41_estimate_noise.sh $SEQUENCE $F $L
cp output_data/1_preprocessing/$SEQUENCE/sigma.txt $INPUT_DIR/

# call different denoisers

# nldct
./41_run_nldct_denoising.sh $SEQUENCE $F $L

# backward nl-kalman
./41_run_nlkalmanbwd_denoising.sh $SEQUENCE $F $L

# recursive bilateral filter
./41_run_rbilf_denoising.sh $SEQUENCE $F $L
