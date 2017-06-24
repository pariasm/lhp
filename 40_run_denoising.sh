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
# ./41_run_nldct_denoising.sh $SEQUENCE $F $L

# nlbayes
./41_run_nlbayes_denoising.sh $SEQUENCE $F $L

# r-nlbayes

# etc...
