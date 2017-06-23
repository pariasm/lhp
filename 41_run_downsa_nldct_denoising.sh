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
echo "Denoising sequence $SEQUENCE. Output stored output_data/4_denoising/$SEQUENCE/nldct/"

INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
OFLOW_DIR="output_data/3_oflow/$SEQUENCE"
OUTPUT_DIR="output_data/4_denoising/$SEQUENCE/nldct/downsa"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# create directory
mkdir -p $OUTPUT_DIR

# downsample the sequence
PLAMBDA="src/utils/imscript/bin/plambda"
DOWNSA="src/utils/imscript/bin/downsa"
ZF=2
for i in $(seq -f "%03g" $F $L)
do
	$DOWNSA v $ZF ${INPUT_DIR}/$i.tif   $OUTPUT_DIR/$i.tif
	$DOWNSA v $ZF ${OFLOW_DIR}/$i.f.flo $OUTPUT_DIR/$i.f.flo
	$DOWNSA v $ZF ${OFLOW_DIR}/$i.b.flo $OUTPUT_DIR/$i.b.flo
	$PLAMBDA $OUTPUT_DIR/$i.f.flo "x $ZF /" -o $OUTPUT_DIR/$i.f.flo
	$PLAMBDA $OUTPUT_DIR/$i.b.flo "x $ZF /" -o $OUTPUT_DIR/$i.b.flo
done

# nldct binary
DENO="src/4_denoising/nldct/build/bin/vnlbayes"
SIGMA=$(cat "output_data/1_preprocessing/$SEQUENCE/sigma.txt")
SIGMA=$(bc <<< "$SIGMA / $ZF")

# run denoising (first step)
$DENO \
-i ${OUTPUT_DIR}/%03d.tif -f $F -l $L -sigma $SIGMA -has-noise \
-fof ${OUTPUT_DIR}/%03d.f.flo -bof ${OUTPUT_DIR}/%03d.b.flo \
-px2 0 -px1 16 -pt1 3 -wx1 31 -wt1 6 -np1 240 -b1 2.0 \
-bsic ${OUTPUT_DIR}/b_%03d.tif

mv measures.txt ${OUTPUT_DIR}/measures_basic

# run denoising (second step)
$DENO \
-i ${OUTPUT_DIR}/%03d.tif -f $F -l $L -sigma $SIGMA -has-noise \
-b ${OUTPUT_DIR}/b_%03d.tif \
-fof ${OUTPUT_DIR}/%03d.f.flo -bof ${OUTPUT_DIR}/%03d.b.flo \
-px1 0 -px2 16 -pt2 3 -wx2 31 -wt2 6 -np2 120 \
-deno ${OUTPUT_DIR}/d_%03d.tif

# clean
rm {bsic,diff}_???.png
rm ${OUTPUT_DIR}/*.flo
mv measures.txt ${OUTPUT_DIR}/measures


