#!/bin/bash
# Evals vnlm using ground truth

SEQ=$1 # sequence path
FFR=$2 # first frame
LFR=$3 # last frame
SIG=$4 # noise standard dev.
OUT=$5 # output folder
PRM=$6 # denoiser parameters

mkdir -p $OUT/s$SIG
OUT=$OUT/s$SIG
 
# add noise {{{1
for i in $(seq $FFR $LFR);
do
	file=$(printf $OUT/"%03d.tif" $i)
	if [ ! -f $file ]
	then
		export SRAND=$RANDOM;
		awgn $SIG $(printf $SEQ $i) $file
	fi
done

# compute optical flow {{{1
TVL1="/home/pariasm/Work/optical_flow/algos/tvl1flow_3/tvl1flow"
for i in $(seq $((FFR+1)) $LFR);
do
	file=$(printf $OUT"/%03d_b.flo" $i)
	if [ ! -f $file ]
	then
		$TVL1 $(printf $OUT"/%03d.tif" $i) \
				$(printf $OUT"/%03d.tif" $((i-1))) \
				$file \
				0 0.25 0.2 0.3 100 0.5 5 0.01 0; 
	fi
done
cp $(printf $OUT"/%03d_b.flo" $((FFR+1))) $(printf $OUT"/%03d_b.flo" $FFR)

# run denoising {{{1
../build/bin/vnlmeans \
 -i $OUT"/%03d.tif" -o $OUT"/%03d_b.flo" -f $FFR -l $LFR -s $SIG \
 -d $OUT"/deno_%03d.tif" $PRM

# compute psnr {{{1
SS=0
for i in $(seq $FFR $LFR);
do
	MM[$i]=$(psnr.sh $(printf $SEQ $i) $(printf $OUT/"deno_%03d.tif" $i) m)
	n=$((i - FFR))
	SS=$(echo "($SS*$n + ${MM[$i]})/$((n+1))" | bc -l)
	MM[$i]=$(plambda -c "${MM[$i]} sqrt")
	PP[$i]=$(plambda -c "255 ${MM[$i]} / log10 20 *")
done

echo "Frame RMSE " ${MM[*]}  > $OUT/measures
echo "Frame PSNR " ${PP[*]} >> $OUT/measures

RMSE=$(plambda -c "$SS sqrt")
PSNR=$(plambda -c "255 $RMSE / log10 20 *")
echo "Total RMSE $RMSE" >> $OUT/measures
echo "Total PSNR $PSNR" >> $OUT/measures


# vim:set foldmethod=marker:
