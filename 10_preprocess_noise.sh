#! /bin/bash

# command line inputs
SEQUENCE=$1 # sequence folder
F=${2:-1}   # first frame (optional: default 1)
L=${3:-0}   # last frame  (optional: default all frames)

# orientation of bands (vertical or horizontal)
BANDS_DIR="vertical"

# number of parallel threads
NUM_PROCS=30

# downsample the sequence as part of the preprocessing
DOWNSAMPLE=0
# -----------------------------------------------------------------------
# NOTE: currently we downsample the sequence AFTER the stabilization
# step. The reason is the following: the stabilized sequence is obtained
# by resampling the original sequence. This resampling uses bicubic
# interpolation, and correlates the noise. We apply the downsampling
# after the stabilization as a way of reducing this correlation.
#
# The drawback is that the stabilization step takes more time. For a
# faster pipeline, set DOWNSAMPLE=1 in this script and set DOWNSAMPLE=0
# after the stabilization (see 20_stabilize_video.sh).
# -----------------------------------------------------------------------

# check correct number of input args
if [ "$#" -lt 1 ]; 
then
	echo "Usage: $0 sequence-folder [first-frame last-frame]" >&2
	exit 1
fi

echo "Preprocessing sequence $SEQUENCE. Output stored in output_data/1_preprocessing/$SEQUENCE/"

# determine last frame
if [ $L -lt 1 ];
then
	L=$(ls input_data/$SEQUENCE/???.tif | wc -l)
fi

# create sym links in the output folder
OUTPUT_DIR="output_data/1_preprocessing/$SEQUENCE"
mkdir -p $OUTPUT_DIR
BASE_DIR=$(pwd)
cd $OUTPUT_DIR
for i in $(seq -f "%03g" $F $L)
do
	ln -s ${BASE_DIR}/input_data/${SEQUENCE}/$i.tif i$i.tif
done
cd -

# add bin folder to path
export PATH=`pwd`/bin/:$PATH

# the preprocessing is written in a makefile
# running make -j takes care of the parallelization
cp preprocess.mk $OUTPUT_DIR/Makefile
make BANDS_DIRECTION=$BANDS_DIR -C $OUTPUT_DIR -j $NUM_PROCS

# change name of the output files to match next scripts
cd $OUTPUT_DIR
for i in $(seq -f "%03g" $F $L)
do
	ln -s s_m_ub_lin_i$i.tif $i.tif
done
cd -

# downsample the sequence
if [ $DOWNSAMPLE -eq 1 ]; 
then
	DOWNSA="src/utils/imscript/bin/downsa"
	ZF=2
	for i in $(seq -f "%03g" $F $L)
	do
		$DOWNSA v $ZF $OUTPUT_DIR/$i.tif $OUTPUT_DIR/$i.tif
	done | parallel
fi

# scale outputs to approx 0 255 range (required for stabilization)

# first compute the 1%-99% range for a subset of the images
S=10 # step (only 1 out of $S images)
for i in $(seq -f "%03g" $F $S $L)
do
	imprintf "%q[1] %q[99]\n" $OUTPUT_DIR/$i.tif >> $OUTPUT_DIR/s_m_ub_lin_iranges.txt
done

# average the previous ranges, to have a unique 
# approximate range for the sequence
awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%f ", a[i]/NR;}; printf "\n"}' \
	$OUTPUT_DIR/s_m_ub_lin_iranges.txt > \
	$OUTPUT_DIR/s_m_ub_lin_irange.txt

# now scale estimated range to 0 255
RANGE=$(cat $OUTPUT_DIR/s_m_ub_lin_irange.txt)
for i in $(seq -f "%03g" $F $L)
do
	plambda $OUTPUT_DIR/$i.tif "$RANGE range 255 *" -o $OUTPUT_DIR/$i.tif
done | parallel

# run ponomarenko's noise sigma estimator
# ---------------------------------------------------------------------
# NOTE: noise estimation should be done before denoising. If the
# downsampling is performed after the stabilization, it would make 
# sense to estimate the noise of the downsampled sequence. However, it
# works better estimating the noise before stabilization. This might be
# because even if the downsampling is applied after the stabilization, the
# noise is still correlated, and this harms the estimation algorithm.
# To speed up the noise estimation, and also to take the future down
# scaling into account, we estimate the noise on downscaled images.
# ---------------------------------------------------------------------
./11_estimate_noise.sh $SEQUENCE $F $L


