#!/bin/bash

#this script is used to execute the program from the IPOL demo
bindir=$1
vidin=$2
denoout=$3
nisyout=$4
sigma=$5
px=$6
wx=$7
hx=$8
hd="$9"
ht=${10}
lambda=${11}
pxlfmt=${12}

#extract info from video
info=`avprobe -v error -show_streams  $vidin`
info="${info#*codec_type=video}"
echo

width=`echo ${info#*width=}| cut -d' ' -f 1` 
height=`echo ${info#*height=}| cut -d' ' -f 1` 
framerate=`echo ${info#*avg_frame_rate=}| cut -d' ' -f 1 | sed 's/\/1//'`
nframes=`echo ${info#*nb_frames=}| cut -d' ' -f 1`
size=${width}x${height}

if [ "$framerate" == "0/0" ] ; then
  echo "Error reading the frame rate of the video (default 30)"
  framerate="30"
fi

echo input parameters and video info
echo " "sigma: $sigma
echo " "psz: $px
echo " "wsz: $wx
echo " "hx : $hx
echo " "hd : $hd
echo " "ht : $hx
echo " "lambda: $lambda
echo " "input video size: $width x $height x $nframes @ $framerate fps
echo

echo extract frames from input video
if [ $pxlfmt == "gray" ]
then
	echo avconv -v error -i $vidin -f image2 -vf format=gray i%04d.png
	time avconv -v error -i $vidin -f image2 -vf format=gray i%04d.png

	for i in $(seq 1 $nframes)
	do
		plambda i$(printf %04d $i).png "x[0]" -o i$(printf %04d $i).png
	done
else
	echo avconv -v error -i $vidin -f image2 i%04d.png
	time avconv -v error -i $vidin -f image2 i%04d.png
fi
echo

echo run rnlmeans denoising
export OMP_NUM_THREADS=24 # set max number of threads
echo $bindir/rbilf-gt.sh i%04d.png 1 $nframes $sigma . \
	  "-p $px -w $wx --whx $hx --whd $hd --wht $ht --lambda $lambda -v"
time $bindir/rbilf-gt.sh i%04d.png 1 $nframes $sigma . \
	  "-p $px -w $wx --whx $hx --whd $hd --wht $ht --lambda $lambda -v"
echo

echo convert tifs to png
for ((i=1; i<=$nframes; i++))
do
	plambda n$(printf "%04d" $i).tif "x" -o n$(printf "%04d" $i).png
	plambda d$(printf "%04d" $i).tif "x" -o d$(printf "%04d" $i).png
done

echo save output video as lossless mp4
echo avconv -y -v error -framerate $framerate -f image2 -i d%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $denoout
time avconv -y -v error -framerate $framerate -f image2 -i d%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $denoout

echo avconv -y -v error -framerate $framerate -f image2 -i n%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $nisyout
time avconv -y -v error -framerate $framerate -f image2 -i n%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $nisyout

if [ $pxlfmt == "gray" ]
then
	# overwrite input video (if it was RGB, it will be overwritten by
	# a grayscale video)
	echo avconv -y -v error -framerate $framerate -f image2 -i i%04d.png \
	     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $vidin
	time avconv -y -v error -framerate $framerate -f image2 -i i%04d.png \
	     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $vidin
fi
echo

nkeep=19
echo keep $nkeep frames at the beginning of the sequence
f1=1
f2=$((nkeep > nframes ? nframes : nkeep))

echo remove tifs, flos and unnecessary pngs
rm -R *.tif *.flo
for ((i=$nkeep+1; i<=$nframes; i++))
do
	rm -R i$(printf "%04d" $i).png
done
