#! /bin/bash

# command line inputs
SEQUENCE=$1 # sequence folder
F=${2:-1}   # first frame (optional: default 1)
L=${3:-0}   # last frame  (optional: default all frames)

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

INPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
SIGMAS="output_data/1_preprocessing/$SEQUENCE/sigmas.txt"
OUTPUT="output_data/1_preprocessing/$SEQUENCE/sigma.txt"

echo "	Estimating noise for sequence $SEQUENCE. Output stored in $OUTPUT"


# determine extension
FF=$(printf %03d $F)
if [ -f "input_data/${SEQUENCE}/${FF}.tif" ]
then
	EXT="tif"
elif [ -f "input_data/${SEQUENCE}/${FF}.tiff" ]
then
	EXT="tiff"
elif [ -f "input_data/${SEQUENCE}/${FF}.png" ]
then
	EXT="png"
else
	echo "File ${SEQUENCE}/${FF}.{tif,tiff,png} not found"
	exit 1
fi

# determine last frame
if [ $L -lt 1 ];
then
	L=$(ls input_data/$SEQUENCE/???.$EXT | wc -l)
fi


PONO=src/1_preprocessing/ponomarenko/ponomarenko
CONVICON=src/utils/convicon/bin/convicon
DOWNSA="src/utils/imscript/bin/downsa"

# downsampling factor
ZF=2

# step (only run ponomarenko in 1/S of the frames)
S=10

for i in $(seq -f "%03g" $F $S $L)
do
	# downsample to simulate the coarse-scale noise
 	$DOWNSA v $ZF $INPUT_DIR/$i.tif $INPUT_DIR/tmp.tif

	# convert from tif to RGB
	$CONVICON -i $INPUT_DIR/tmp.tif -o $INPUT_DIR/tmp.RGB

	# run ponomarenko's noise estimator with a single bin
	$PONO -b 1 $INPUT_DIR/tmp.RGB | awk '{print $2}' >> $SIGMAS
done | parallel

# compute average sigma
awk '{s+=$1} END {print s/NR}' RS=" " $SIGMAS > $OUTPUT
rm -f $INPUT_DIR/tmp.RGB $INPUT_DIR/tmp.RGB

