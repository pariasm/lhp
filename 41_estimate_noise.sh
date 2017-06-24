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
INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
SIGMAS="$INPUT_DIR/sigmas.txt"
OUTPUT="$INPUT_DIR/sigma.txt"

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
	L=$(ls $INPUT_DIR/???.$EXT | wc -l)
fi


PONO=src/1_preprocessing/ponomarenko/ponomarenko
CONVICON=src/utils/convicon/bin/convicon

# step (only run ponomarenko in 1/S of the frames)
S=10

for i in $(seq -f "%03g" $F $S $L)
do
	# convert from tif to RGB
	$CONVICON -i $INPUT_DIR/$i.tif -o $INPUT_DIR/tmp.RGB

	# run ponomarenko's noise estimator with a single bin
	$PONO -b 1 $INPUT_DIR/tmp.RGB | awk '{print $2}' >> $SIGMAS
done | parallel

# compute average sigma
awk '{s+=$1} END {print s/NR}' RS=" " $SIGMAS > $OUTPUT
rm -f $INPUT_DIR/tmp.RGB

