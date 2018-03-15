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
vnlmeans
vnlm-gt.sh
vnlm-mp4.sh
psnr.sh
```

vnlmeans is the denoising algorithm. Run `vnlmeans -h` for more help. The
scripts are necessary because vnlmeans only does denoising, meanining it
doesn't compute the optical flow (which optional but highly recommended), add
noise or compute the psnr. 

`vnlm-mp4.sh`:	call vnlm-gt.sh but receiving and output mp4 videos. Requires
`avconv` or `ffmpeg` to be installed.

`vnlm-gt.sh`: run vnlmeans and compare the output with the ground truth.
1. adds noise (requires `awgn` from imscript)
2. computes optical flow (requires `downsa` and `upsa` from impscript)
3. computes the psnr of the output (requires from impscript `plambda`, `crop`,
	`imprintf`)
4. leaves all data in an output folder which can be given as a parameter (noisy
	frames, optical flow, denoised, and quantitative measurements in `measures.txt`)

Example
-------

```
$ bin/vnlm-gt.sh \
		/path/to/sequence/%03d.png \
		first-frame last-frame sigma output-dir \
		"string with the params for vnlmeans"
```

For example:
```
bin/vnlm-gt.sh ~/my-seq/%03d.png 1 30 20 parkjoy2 "-p 8 -w 10 -v"
```


This will (first add noise) and denoise a sequence of 30 png files.
It will produce the following output to stdout:

```
/home/user/projects/rnlmeans/build/bin/vnlmeans -i parkjoy2/n%04d.tif -o parkjoy2/%04d_b.flo -f 1 -l 30 -s 20 -d parkjoy2/d%04d.tif -p 8 -w 10 -v
parameters:
    noise  20.000000
    patch-wise mode
    patch     8
    search    10
    w_hx      17
    w_hd      -nan
    w_thx     0.05
    w_ht      27
    w_htv     27
    lambda    1
    tv_lambda 0.85

WEIGHTED_AGGREGATION ON
AGGREGATE_TRANSITION_VAR ON
loading video parkjoy2/n%04d.tif
loading flow parkjoy2/%04d_b.flo
processing frame 1
processing frame 2
processing frame 3
...
(etcera)
```
