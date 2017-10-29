#!/bin/bash
# Tune the algorithm's parameters

# noise levels
sigmas=(10 20 40 30)

# fixed parameters
pszs=(4 8 12)
wszs=(5 10 15)

# number of trials
ntrials=100

# test sequences
seqs=(\
derf/bus_mono \
derf/foreman_mono \
derf/football_mono \
derf/tennis_mono \
derf/stefan_mono \
tut/gsalesman \
)

# seq folder
sf='/home/pariasm/Remote/avocat/denoising/data/'

for ((i=0; i < $ntrials; i++))
do
	# randomly draw noise level and parameters

	r=$(awk -v M=2 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
	s=${sigmas[$r]}

#	r=$(awk -v M=2 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
#	p=${pszs[$r]}

#	r=$(awk -v M=2 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
#	w=${wszs[$r]}

	p=8
	w=10

	W=$((4*s*s))
	whx=$(awk -v M=$W -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
	wht=$(awk -v M=$W -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
#	whv=$(awk -v M=$W -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M+1))}')
	whv=0

	lambda=$(awk -v s=$RANDOM 'BEGIN{srand(s); print rand()}')

	trialfolder=$(printf "trials/rnlm.s%02d.p%02d.w%02d.whx%04d.wht%04d.whv%04d.lambda%5.3f\n" \
		$s $p $w $whx $wht $whv $lambda)

	params=$(printf " -p %d -w %d --whx %d --wht %d --whtv %d --lambda %f" \
		$p $w $whx $wht $whv $lambda)

	mpsnr=0
	nseqs=${#seqs[@]}
	nf=10
	if [ ! -d $trialfolder ]
	then
		for seq in ${seqs[@]}
		do
			echo "./vnlm_train.sh ${sf}${seq} 1 $nf $s $trialfolder \"$params\""
			./vnlm_train.sh ${sf}${seq} 1 $nf $s $trialfolder "$params"
			psnr=$(./vnlm_train.sh ${sf}${seq} 1 $nf $s $trialfolder "$params")
			mpsnr=$(echo "$mpsnr + $psnr/$nseqs" | bc -l)
			#echo $mpsnr
		done
	fi
	
	printf "%2d %2d %2d %4d %4d %4d %5.3f %7.4f\n" \
		$s $p $w $whx $wht $whv $lambda $mpsnr >> trials/table

	rm $trialfolder/*.tif

done

#	r=$(awk -v m=1 -v M=10 -v s=$RANDOM 'BEGIN{srand(s); print int(m+rand()*(M-m+1))}')
