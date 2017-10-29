#!/bin/bash
# This script is used by train.sh. It uses a set of precomputed 
# noisy sequences with their optical flow. To speed up training,
# the noise is already added, and the optical flow is already 
# computed.

SEQ=$1 # sequence path
FFR=$2 # first frame
LFR=$3 # last frame
SIG=$4 # noise standard dev.
OUT=$5 # output folder
PRM=$6 # denoiser parameters

# output folder
mkdir -p $OUT
 
# folders with pre-computed data
ORIG="$SEQ/%03d.png"
NISY="$SEQ/s${SIG}/%03d.tif"
FLOW="$SEQ/s${SIG}/tvl1_%03d_b.flo"

# run denoising
../build/bin/vnlmeans \
 -i $NISY -o $FLOW -f $FFR -l $LFR -s $SIG \
 -d $OUT"/deno_%03d.tif" $PRM

# compute psnr
SS=0
for i in $(seq $FFR $LFR);
do
	MM[$i]=$(psnr.sh $(printf $ORIG $i) $(printf $OUT/"deno_%03d.tif" $i) m)
	n=$((i - FFR))
	SS=$(echo "($SS*$n + ${MM[$i]})/$((n+1))" | bc -l)
	MM[$i]=$(plambda -c "${MM[$i]} sqrt")
	PP[$i]=$(plambda -c "255 ${MM[$i]} / log10 20 *")
done
 
echo "Frame RMSE " ${MM[*]}  > $OUT/$(basename $SEQ).measures
echo "Frame PSNR " ${PP[*]} >> $OUT/$(basename $SEQ).measures
 
RMSE=$(plambda -c "$SS sqrt")
PSNR=$(plambda -c "255 $RMSE / log10 20 *")
echo "Total RMSE $RMSE" >> $OUT/$(basename $SEQ).measures
echo "Total PSNR $PSNR" >> $OUT/$(basename $SEQ).measures
echo $PSNR


# vim:set foldmethod=marker:
