#! /bin/bash

BANDS_DIR="vertical"

NUM_PROCS=30
export PATH=`pwd`/bin/:$PATH

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

F=${2:-1} # first frame
L=${3:-0} # last frame

SEQUENCE=$1
echo "Preprocessing sequence $SEQUENCE. Output stored in output_data/1_preprocessing/$SEQUENCE/"

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


# for now, we just create sym links in the output folder
OUTPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
mkdir -p $OUTPUT_DIR
BASE_DIR=$(pwd)
cd $OUTPUT_DIR
for i in $(seq -f "%03g" $F $L)
do
	ln -s ${BASE_DIR}/input_data/${SEQUENCE}/$i.$EXT i$i.$EXT
done
cd -

cp preprocess.mk $OUTPUT_DIR/Makefile
make BANDS_DIRECTION=$BANDS_DIR -C $OUTPUT_DIR -j $NUM_PROCS

# change name of enric's output files to match next scripts
cd $OUTPUT_DIR
for i in $(seq -f "%03g" $F $L)
do
	ln -s s_m_ub_lin_i$i.$EXT $i.$EXT
done
cd -

# # downsample the sequence
# DOWNSA="src/utils/imscript/bin/downsa"
# ZF=2
# for i in $(seq -f "%03g" $F $L)
# do
# 	$DOWNSA v $ZF $OUTPUT_DIR/$i.tif $OUTPUT_DIR/$i.tif
# done | parallel

# scale outputs to approx 0 255 range (required for stabilization)

# first compute the 1%-99% range for a bunch of images
S=10 # step
for i in $(seq -f "%03g" $F $S $L)
do
	imprintf "%q[1] %q[99]\n" $OUTPUT_DIR/$i.$EXT >> $OUTPUT_DIR/s_m_ub_lin_iranges.txt
done

awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%f ", a[i]/NR;}; printf "\n"}' \
	$OUTPUT_DIR/s_m_ub_lin_iranges.txt > \
	$OUTPUT_DIR/s_m_ub_lin_irange.txt


RANGE=$(cat $OUTPUT_DIR/s_m_ub_lin_irange.txt)
for i in $(seq -f "%03g" $F $L)
do
	plambda $OUTPUT_DIR/$i.$EXT "$RANGE range 255 *" -o $OUTPUT_DIR/$i.$EXT
done | parallel

# run ponomarenko's noise estimator
./11_estimate_noise.sh $SEQUENCE $F $L


