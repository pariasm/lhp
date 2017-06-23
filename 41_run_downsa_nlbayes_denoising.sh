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
echo "Denoising sequence $SEQUENCE. Output stored output_data/4_denoising/$SEQUENCE/nlbayes/"

INPUT_DIR="output_data/2_stabilization/$SEQUENCE"
OFLOW_DIR="output_data/3_oflow/$SEQUENCE"
OUTPUT_DIR="output_data/4_denoising/$SEQUENCE/nlbayes/downsa"

# determine last frame
if [ $L -lt 1 ];
then
	N=$(ls $INPUT_DIR/???.tif | wc -l)
	L=$((F + N - 1))
fi

# nlbayes binary
DENO="src/4_denoising/nlbayes/build/bin/vnlbayes"
SIGMA=$(cat "output_data/1_preprocessing/$SEQUENCE/sigma.txt")
ZF=2
SIGMA=$(bc -l <<< ${SIGMA}/${ZF})

echo $SIGMA

# create directory
mkdir -p $OUTPUT_DIR

# downsample the sequence
DOWNSA="src/utils/imscript/bin/downsa"
PLAMBDA="src/utils/imscript/bin/plambda"
for i in $(seq -f "%03g" $F $L)
do
	$DOWNSA v $ZF ${INPUT_DIR}/$i.tif $OUTPUT_DIR/$i.tif
	$DOWNSA v $ZF ${OFLOW_DIR}/$i.f.flo $OUTPUT_DIR/$i.f.flo
	$DOWNSA v $ZF ${OFLOW_DIR}/$i.b.flo $OUTPUT_DIR/$i.b.flo
	$PLAMBDA $OUTPUT_DIR/$i.f.flo x\ $ZF\ / -o $OUTPUT_DIR/$i.f.flo
	$PLAMBDA $OUTPUT_DIR/$i.b.flo x\ $ZF\ / -o $OUTPUT_DIR/$i.b.flo
done

# run denoising (first step)
$DENO \
-i ${OUTPUT_DIR}/%03d.tif -f $F -l $L -sigma $SIGMA -has-noise \
-fof $OUTPUT_DIR/%03d.f.flo -bof $OUTPUT_DIR/%03d.b.flo \
-px2 0 -px1 12 -pt1 3 -wx1 31 -wt1 6 -r1 60 -np1 200 -th1 2.1 \
-bsic ${OUTPUT_DIR}/b_%03d.tif

mv measures.txt ${OUTPUT_DIR}/measures_basic

# run denoising (second step)
$DENO \
-i ${OUTPUT_DIR}/%03d.tif -f $F -l $L -sigma $SIGMA -has-noise \
-b ${OUTPUT_DIR}/b_%03d.tif \
-fof $OUTPUT_DIR/%03d.f.flo -bof $OUTPUT_DIR/%03d.b.flo \
-px1 0 -px2 12 -pt2 3 -wx2 31 -wt2 6 -r2 60 -np2 100 -th2 1.45 -flat-area2 \
-deno ${OUTPUT_DIR}/d_%03d.tif

# clean
rm {bsic,diff}_???.png
mv measures.txt ${OUTPUT_DIR}/measures


