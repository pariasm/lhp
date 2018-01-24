#! /bin/bash

# check correct number of input args
if [ "$#" -ne 4 ]; 
then
	echo "Usage: $0 sequence-path first-frame last-frame output-folder" >&2
	exit 1
fi

# Computes tvl1 optical flow for a (noisy) sequence. 
D=${1:-"example/sequence/path/frame_%03d.tif"}  # sequence path, in printf format (MUST have a %d identifier)
F=${2:-1}                                       # first frame
L=${3:-150}                                     # last frame 
O=${4:-"."}                                     # output folder

echo "	Computing TVL1 flow for sequence $D. Output stored in $O"


TVL1="src/3_oflow/tvl1flow_3/tvl1flow"

mkdir -p $O

DW=0.1 # weight of data attachment term 0.1 ~ very smooth 0.2 ~ noisy
FS=1   # finest scale (0 image scale, 1 one coarser level, 2 more coarse, etc...

# compute forward flow
for i in `seq $F $((L - 1))`;
do
	$TVL1 `printf ${D} $i` `printf ${D} $((i + 1))` \
		   ${O}/`printf %03d.f.flo $i` \
		   0 0.25 $DW 0.3 100 $FS 0.5 5 0.01 0; 
done
cp ${O}/`printf %03d.f.flo $((L - 1))` ${O}/`printf %03d.f.flo $L`

# compute backward flow
for i in `seq $((F + 1)) $L`;
do
	$TVL1 `printf ${D} $i` `printf ${D} $((i - 1))` \
		   ${O}/`printf %03d.b.flo $i` \
		   0 0.25 $DW 0.3 100 $FS 0.5 5 0.01 0; 
done
cp ${O}/`printf %03d.b.flo $((F + 1))` ${O}/`printf %03d.b.flo $F`

