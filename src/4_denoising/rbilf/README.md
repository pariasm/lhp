Video denoising with a recursive bilateral filter
=================================================

Implementation of a recursive version of the bilateral filter for video denoising.

* Author    : Pablo Arias <pariasm@gmail.com>
* Copyright : (C) 2018 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt


Building
--------

Building has been tested on Ubuntu Unix 16.04.

Dependencies:
 * libtiff, libpng, libjpeg

Apart from these, some scripts provided by the code use some functions of 
[imscript](https://github.com/mnhrdt/imscript). These functions need to be
accesible from the shell path.

This code ships with the following 3rd party libs:
 * [iio](https://github.com/mnhrdt/iio)
 * [tvl1](https://doi.org/10.5201/ipol.2013.26)
 * [argparse](https://github.com/cofyc/argparse)

 
Configure and compile the source code using cmake and make. 
For example:

```
$ mkdir build; cd build
$ cmake ..
$ make
```

Usage
--------

After building you'll find in the build/bin/ the following files:

```
tvl1flow
rbilf
rbilf-gt.sh
rbilf-mp4.sh
psnr.sh
```

rbilf is the denoising algorithm. Run `rbilf -h` for more help. The
scripts are necessary because rbilf only does denoising, meanining it
doesn't compute the optical flow (which optional but highly recommended), add
noise or compute the psnr. 

`rbilf-mp4.sh`:	call rbilf-gt.sh but receiving and output mp4 videos. Requires
`avconv` or `ffmpeg` to be installed.

`rbilf-gt.sh`: run rbilf and compare the output with the ground truth.
1. adds noise (requires `awgn` from imscript)
2. computes optical flow (requires `downsa` and `upsa` from impscript)
3. computes the psnr of the output (requires from impscript `plambda`, `crop`,
	`imprintf`)
4. leaves all data in an output folder which can be given as a parameter (noisy
	frames, optical flow, denoised, and quantitative measurements in `measures.txt`)

Example
-------

```
$ bin/rbilf-gt.sh \
		/path/to/sequence/%03d.png \
		first-frame last-frame sigma output-dir \
		"string with the params for rbilf"
```

For example:
```
bin/rbilf-gt.sh /path/to/sequence/bus/%03d.png 1 15 20 bus "-w 3 -v --whx 30"
```

This will add noise to the sequence bus, compute the optical flow and then
denoise it. The output will be stored in the bus folder.
It will produce the following output to stdout:

```
/home/pariasm/Work/denoising/projects/rbilf/build/bin/rbilf -i bus/n%04d.tif -o bus/%04d_b.flo -f 1 -l 15 -s 20 -d bus/d%04d.tif -v --whx 30 -w 3
parameters:
	noise  20.000000
	search    3
	w_hx      30
	w_hd      1.6
	w_thx     0.05
	w_ht      14
	lambda_x  0.1
	lambda_t  0.5

loading video bus/n%04d.tif
loading flow bus/%04d_b.flo
processing frame 1
processing frame 2
processing frame 3
processing frame 4
processing frame 5
processing frame 6
processing frame 7
processing frame 8
processing frame 9
processing frame 10
processing frame 11
processing frame 12
processing frame 13
processing frame 14
processing frame 15

```
