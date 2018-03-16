#!/bin/bash
# Tune the algorithm's parameters

# command to launch in parallel
# $ for m in vari novari; do for i in {00..13}; do echo "./train-$m.sh $m$i 2>/dev/null"; done; done | parallel

# noise levels
sigmas=(10 20 40 30)

# fixed parameters
pszs=(4 8 12)
wszs=(5 10 15)

# number of trials
ntrials=2000

# test sequences
# seqs=(\
# derf/bus_mono \
# derf/foreman_mono \
# derf/football_mono \
# derf/tennis_mono \
# derf/stefan_mono \
# )
# #tut/gsalesman \
seqs=(\
derf-hd/park_joy \
derf-hd/speed_bag \
derf-hd/station2 \
derf-hd/sunflower \
derf-hd/tractor \
)

# seq folder
sf='/home/pariasm/denoising/data/'

output=${1:-"trials"}

for ((i=0; i < $ntrials; i++))
do
	# randomly draw noise level and parameters

	# noise level
	r=$(awk -v M=2 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
	s=${sigmas[$r]}

	# patch size
#	r=$(awk -v M=2 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
#	p=${pszs[$r]}

	# search region
#	r=$(awk -v M=2 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
#	w=${wszs[$r]}

	p=4
	w=10

	# spatial and temporal weights
	W=$((2*s))
	whx=$(awk -v M=$W -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
	wht=$(awk -v M=$W -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
	whv=$(awk -v M=$W -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
#	whv=0

	lambda=$(awk -v s=$RANDOM 'BEGIN{srand(s); print rand()}')
	lambtv=$(awk -v s=$RANDOM 'BEGIN{srand(s); print rand()}')

	trialfolder=$(printf "$output/rnlm.s%02d.p%02d.w%02d.whx%04d.wht%04d.whv%04d.lambtv%5.3f.lambda%5.3f\n" \
		$s $p $w $whx $wht $whv $lambtv $lambda)

	params=$(printf " -p %d -w %d --whx %d --wht %d --whtv %d --lambtv %f --lambda %f" \
		$p $w $whx $wht $whv $lambtv $lambda)

	echo $trialfolder

	mpsnr=0
	nseqs=${#seqs[@]}
	f0=70
	f1=85
	if [ ! -d $trialfolder ]
	then
		for seq in ${seqs[@]}
		do
			echo "./rbilf_train.sh ${sf}${seq} $f0 $f1 $s $trialfolder \"$params\""
			psnr=$(./rbilf_train.sh ${sf}${seq} $f0 $f1 $s $trialfolder "$params")
			mpsnr=$(echo "$mpsnr + $psnr/$nseqs" | bc -l)
			#echo $mpsnr
		done
	fi
	
	printf "%2d %2d %2d %4d %4d %4d %5.3f %5.3f %7.4f\n" \
		$s $p $w $whx $wht $whv $lambtv $lambda $mpsnr >> $output/table

	rm $trialfolder/*.tif

done

#	r=$(awk -v m=1 -v M=10 -v s=$RANDOM 'BEGIN{srand(s); print int(m+rand()*(M-m+1))}')
