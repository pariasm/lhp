#include "argparse.h"   // command line parser
#include "iio.h"        // image i/o

#include <assert.h>     // assert
#include <stdlib.h>     // NULL, EXIT_SUCCESS/FAILURE
#include <stdio.h>      // printf, fprintf, stderr
#include <math.h>       // expf, nans (used as boundary value by bicubic interp)
#include <string.h>

// some macros and data types [[[1

#define FLT_HUGE 1e10

// read/write image sequence [[[1
static
float * vio_read_video_float_vec(const char * const path, int first, int last, 
		int *w, int *h, int *pd)
{
	// retrieve size from first frame and allocate memory for video
	int frames = last - first + 1;
	int whc;
	float *vid = NULL;
	{
		char frame_name[512];
		sprintf(frame_name, path, first);
		float *im = iio_read_image_float_vec(frame_name, w, h, pd);

		// size of a frame
		whc = *w**h**pd;

		vid = malloc(frames*whc*sizeof(float));
		memcpy(vid, im, whc*sizeof(float));
		if(im) free(im);
	}

	// load video
	for (int f = first + 1; f <= last; ++f)
	{
		int w1, h1, c1;
		char frame_name[512];
		sprintf(frame_name, path, f);
		float *im = iio_read_image_float_vec(frame_name, &w1, &h1, &c1);

		// check size
		if (whc != w1*h1*c1)
		{
			fprintf(stderr, "Size missmatch when reading frame %d\n", f);
			if (im)  free(im);
			if (vid) free(vid);
			return NULL;
		}

		// copy to video buffer
		memcpy(vid + (f - first)*whc, im, whc*sizeof(float));
		if(im) free(im);
	}

	return vid;
}

static
void vio_save_video_float_vec(const char * const path, float * vid, 
		int first, int last, int w, int h, int c)
{
	const int whc = w*h*c;
	for (int f = first; f <= last; ++f)
	{
		char frame_name[512];
		sprintf(frame_name, path, f);
		float * im = vid + (f - first)*whc;
		iio_save_image_float_vec(frame_name, im, w, h, c);
	}
}


// bicubic interpolation [[[1

#ifdef NAN
// extrapolate by nan
inline static
float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
	assert(l >= 0 && l < pd);
	return (i < 0 || i >= w || j < 0 || j >= h) ? NAN : x[(i + j*w)*pd + l];
}
#endif//NAN

inline static
float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static
float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

static
void bicubic_interpolation_nans(float *result,
		float *img, int w, int h, int pd, float x, float y)
{
	x -= 1;
	y -= 1;

	int ix = floor(x);
	int iy = floor(y);
	for (int l = 0; l < pd; l++) {
		float c[4][4];
		for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = getsample_nan(img, w, h, pd, ix + i, iy + j, l);
		float r = bicubic_interpolation_cell(c, x - ix, y - iy);
		result[l] = r;
	}
}

static
void warp_bicubic_inplace(float *imw, float *im, float *of, int w, int h, int ch)
{
	// warp previous frame
	for (int y = 0; y < h; ++y)
	for (int x = 0; x < w; ++x)
	{
		float xw = x + of[(x + y*w)*2 + 0];
		float yw = y + of[(x + y*w)*2 + 1];
		bicubic_interpolation_nans(imw + (x + y*w)*ch, im, w, h, ch, xw, yw);
	}
	return;
}

// recursive nl-means parameters [[[1

// struct for storing the parameters of the algorithm
struct vnlmeans_params
{
	int search_sz;       // search window radius
	float weights_hx;    // spatial patch similarity weights parameter
	float weights_hd;    // spatial distance weights parameter (for bilateral)
	float weights_thx;   // spatial weights threshold
	float weights_ht;    // temporal patch similarity weights parameter
	float dista_lambda;  // weight of current frame in patch distance
	float tv_lambda;     // weight of current frame in patch distance
};

// set default parameters as a function of sigma
void vnlmeans_default_params(struct vnlmeans_params * p, float sigma)
{
	// the parameters are based on the following parameters
	// found with a parameter search:
	//
	//           wsz hx   hd  ht   hv ltv  lx
	// sigma 10:  10 10.2 4    7.5 0  0.5  0.1 
	// sigma 20:  10 24.4 1.6 14.1 0  0.3  0.1
	// sigma 40:  10 48.0 1.6 27.1 0  0.6  0.02
	if (p->weights_hx   < 0) p->weights_hx   = 1.2 * sigma;
	if (p->weights_hd   < 0) p->weights_hd   = 1.6;
	if (p->search_sz    < 0) p->search_sz    = 3*p->weights_hd;
	if (p->weights_thx  < 0) p->weights_thx  = .05f;
	if (p->weights_ht   < 0) p->weights_ht   = 0.7 * sigma;
	if (p->tv_lambda    < 0) p->tv_lambda    = 0.5;
	if (p->dista_lambda < 0)
		p->dista_lambda = fmax(0, fmin(0.2, 0.1 - (sigma - 20)/400));

	// limit search region to prevent too long running times
	p->search_sz = fmin(15, fmin(3*p->weights_hd, p->search_sz));
}

// recursive bilateral filter for frame t [[[1
/* This is a particular case of the vnlmeans_kalman_frame in which we assume
 * that the patch size is 1. In this implementation we removed all loops over
 * the patch domain, and as a result is slightly faster. This version won't be
 * mantained for the moment. The most up-to-date code will be the
 * vnlmeans_kalman_frame. */
void rbilateral_filter_frame(float *deno1, float *nisy1, float *deno0,
		int w, int h, int ch, float sigma,
		const struct vnlmeans_params prms, int frame)
{
	// definitions [[[2

	const float weights_hx2  = prms.weights_hx * prms.weights_hx;
	const float weights_hd2  = prms.weights_hd * prms.weights_hd * 2;
	const float weights_ht2  = prms.weights_ht * prms.weights_ht;
	const float sigma2 = sigma * sigma;

	// set output and aggregation weights to 0
	for (int i = 0; i < w*h*ch; ++i) deno1[i] = 0.;

	// wrap images with nice pointers to vlas
	float (*d1)[w][ch] = (void *)deno1;       // denoised frame t (output)
	const float (*d0)[w][ch] = (void *)deno0; // denoised frame t-1
	const float (*n1)[w][ch] = (void *)nisy1; // noisy frame at t

	// loop on image pixels [[[2
	#pragma omp parallel for
	for (int py = 0; py < h; ++py)
	for (int px = 0; px < w; ++px)
	{
		// spatial average: loop on search region [[[3
		float D1[ch]; // denoised pixel at p in frame t
		for (int c = 0; c < ch ; ++c)	D1[c] = 0;
		float cp = 0.; // sum of similarity weights, used for normalization
		if (weights_hx2)
		{
			const int wsz = prms.search_sz;
			const int wx[2] = {fmax(px - wsz, 0), fmin(px + wsz + 1, w)};
			const int wy[2] = {fmax(py - wsz, 0), fmin(py + wsz + 1, h)};
			for (int qy = wy[0]; qy < wy[1]; ++qy)
			for (int qx = wx[0]; qx < wx[1]; ++qx)
			{
				// compute patch distance [[[4
				float ww = 0;
				const float l = prms.dista_lambda;
				if (d0 && l != 1 && !isnan(d0[qy][qx][0]) && !isnan(d0[py][px][0]))
					// use noisy and denoised patches from previous frame
					for (int c = 0; c < ch ; ++c)
					{
						const float eN1 = n1[qy][qx][c] - n1[py][px][c];
						const float eD0 = d0[qy][qx][c] - d0[py][px][c];
						ww += l * eN1 * eN1 + (1 - l) * eD0 * eD0;
					}
				else
					// use only noisy from previous frame
					for (int c = 0; c < ch ; ++c)
					{
						const float eN1 = n1[qy][qx][c] - n1[py][px][c];
						ww += eN1 * eN1;
					}

				// compute spatial similarity weight ]]]4[[[4
				if (weights_hx2)
					ww = expf(-1 / weights_hx2 * ww / (float)(ch));
				else
					ww = (qx == px && qy == py) ? 1. : 0.;

				if (weights_hd2 < FLT_HUGE)
				{
					if (weights_hd2)
					{
						const float dx = (px - qx), dy = (py - qy);
						ww *= expf(-1 / weights_hd2 * (dx * dx + dy * dy) );
					}
					else
						ww = (qx == px && qy == py) ? 1. : 0.;
				}

				// accumulate on output pixel ]]]4[[[4
				if (ww > prms.weights_thx)
				{
					cp += ww;
					for (int c = 0; c < ch; ++c)
						D1[c] += n1[qy][qx][c] * ww;
				}// ]]]4
			}
		}
		else
		{
			// copy noisy pixel to output
			cp = 1.;
			for (int c = 0; c < ch; ++c)
				D1[c] = n1[py][px][c];
		}

		// store denoised pixel on output image
		float icp = 1. / fmax(cp, 1e-6);
		for (int c = 0; c < ch ; ++c )
			d1[py][px][c] = D1[c] * icp;


		// temporal average with frame t-1 [[[3
		if (d0 && !isnan(d0[py][px][0]))
		{
			// estimate transition variance
			float v01xy = 0.;
			const float l01 = prms.tv_lambda;
			for (int c = 0; c < ch; ++c)
			{
				const float eD = d1[py][px][c] - d0[py][px][c];
				const float eN = n1[py][px][c] - d0[py][px][c];
				v01xy += l01 * fmax(eN * eN - sigma2, 0.f) + (1 - l01) * eD * eD;
			}

			// normalize variance estimate
			v01xy /= (float)ch;

			// compute variance multiplier
			const float w01xy = fmin(1, fmax(0, expf( - 1./weights_ht2 * v01xy )));

//			// no variances: we assume that v0 = v1
//			const float f = fmin(1., fmax(0., 1./(1. + w01xy)));
			// no variances: we assume that v0 = (1 - w01)v1
			const float f = 1. - w01xy;

			// update pixel value
			for (int c = 0; c < ch; ++c)
				d1[py][px][c] = d0[py][px][c] * (1 - f) + d1[py][px][c] * f;
		}
	}

	return; // ]]]2
}


// main funcion [[[1

// 'usage' message in the command line
static const char *const usages[] = {
	"vnlmeans [options] [[--] args]",
	"vnlmeans [options]",
	NULL,
};

int main(int argc, const char *argv[])
{
	// parse command line [[[2

	// command line parameters and their defaults
	const char *nisy_path = NULL;
	const char *deno_path = NULL;
	const char *flow_path = NULL;
	int fframe = 0, lframe = -1;
	float sigma = 0.f;
	bool verbose = false;
	struct vnlmeans_params prms;
	prms.search_sz    = -1;
	prms.weights_hx   = -1.; // -1 means automatic value
	prms.weights_hd   = -1.;
	prms.weights_thx  = -1.;
	prms.weights_ht   = -1.;
	prms.dista_lambda = -1.;
	prms.tv_lambda    = -1.;

	// configure command line parser
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Algorithm options"),
		OPT_STRING ('i', "nisy"  , &nisy_path, "noisy input path (printf format)"),
		OPT_STRING ('o', "flow"  , &flow_path, "backward flow path (printf format)"),
		OPT_STRING ('d', "deno"  , &deno_path, "denoised output path (printf format)"),
		OPT_INTEGER('f', "first" , &fframe, "first frame"),
		OPT_INTEGER('l', "last"  , &lframe , "last frame"),
		OPT_FLOAT  ('s', "sigma" , &sigma, "noise standard dev"),
		OPT_INTEGER('w', "search", &prms.search_sz, "search region radius"),
		OPT_FLOAT  ( 0 , "whx"   , &prms.weights_hx, "spatial patch sim. weights param"),
		OPT_FLOAT  ( 0 , "whd"   , &prms.weights_hd, "spatial distance weights param"),
		OPT_FLOAT  ( 0 , "wthx"  , &prms.weights_thx, "spatial weights threshold"),
		OPT_FLOAT  ( 0 , "wht"   , &prms.weights_ht, "temporal patch sim. weights param"),
		OPT_FLOAT  ( 0 , "lambda", &prms.dista_lambda, "noisy patch weight in patch distance"),
		OPT_FLOAT  ( 0 , "lambtv", &prms.tv_lambda, "noisy patch weight in transition variance"),
		OPT_GROUP("Program options"),
		OPT_BOOLEAN('v', "verbose", &verbose, "verbose output"),
		OPT_END(),
	};

	// parse command line
	struct argparse argparse;
	argparse_init(&argparse, options, usages, 0);
	argparse_describe(&argparse, "\nA video denoiser based on non-local means.", "");
	argc = argparse_parse(&argparse, argc, argv);

	// default value for noise-dependent params
	vnlmeans_default_params(&prms, sigma);

	// print parameters
	if (verbose)
	{
		printf("parameters:\n");
		printf("\tnoise  %f\n", sigma);
		printf("\tsearch    %d\n", prms.search_sz);
		printf("\tw_hx      %g\n", prms.weights_hx);
		printf("\tw_hd      %g\n", prms.weights_hd);
		printf("\tw_thx     %g\n", prms.weights_thx);
		printf("\tw_ht      %g\n", prms.weights_ht);
		printf("\tlambda    %g\n", prms.dista_lambda);
		printf("\ttv_lambda %g\n", prms.tv_lambda);
		printf("\n");
	}

	// load data [[[2
	if (verbose) printf("loading video %s\n", nisy_path);
	int w, h, c; //, frames = lframe - fframe + 1;
	float * nisy = vio_read_video_float_vec(nisy_path, fframe, lframe, &w, &h, &c);
	{
		if (!nisy)
			return EXIT_FAILURE;
	}

	// load optical flow
	float * flow = NULL;
	if (flow_path)
	{
		if (verbose) printf("loading flow %s\n", flow_path);
		int w1, h1, c1;
		flow = vio_read_video_float_vec(flow_path, fframe, lframe, &w1, &h1, &c1);

		if (!flow)
		{
			if (nisy) free(nisy);
			return EXIT_FAILURE;
		}

		if (w*h != w1*h1 || c1 != 2)
		{
			fprintf(stderr, "Video and optical flow size missmatch\n");
			if (nisy) free(nisy);
			if (flow) free(flow);
			return EXIT_FAILURE;
		}
	}

	// run denoiser [[[2
	const int whc = w*h*c, wh2 = w*h*2;
	float * deno = nisy;
	float * warp0 = malloc(whc*sizeof(float));
	float * deno1 = malloc(whc*sizeof(float));
	for (int f = fframe; f <= lframe; ++f)
	{
		if (verbose) printf("processing frame %d\n", f);

		// TODO compute optical flow if absent

		// warp previous denoised frame
		if (f > fframe)
		{
			float * deno0 = deno + (f - 1 - fframe)*whc;
			if (flow)
			{
				float * flow0 = flow + (f - 0 - fframe)*wh2;
				warp_bicubic_inplace(warp0, deno0, flow0, w, h, c);
			}
			else
				// copy without warping
				memcpy(warp0, deno0, whc*sizeof(float));
		}

		// run denoising
		float *nisy1 = nisy + (f - fframe)*whc;
		float *deno0 = (f > fframe) ? warp0 : NULL;
		rbilateral_filter_frame(deno1, nisy1, deno0, w, h, c, sigma, prms, f);

		memcpy(nisy1, deno1, whc*sizeof(float));
	}

	// save output [[[2
	vio_save_video_float_vec(deno_path, deno, fframe, lframe, w, h, c);

	if (deno1) free(deno1);
	if (warp0) free(warp0);
	if (nisy) free(nisy);
	if (flow) free(flow);

	return EXIT_SUCCESS; // ]]]2
}

// vim:set foldmethod=marker:
// vim:set foldmarker=[[[,]]]:
