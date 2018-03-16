// biharmonic solver
// note: timestep must be negative between -0.08 and -0.099

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// utility function that always returns a valid pointer to memory
static void *xmalloc(size_t n)
{
	void *new = malloc(n);
	if (!new)
	{
		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
		exit(1);
	}
	return new;
}

// the type of a "getpixel" function
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by nearest value (useful for Neumann boundary conditions)
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// evaluate the laplacian of image x at point i, j
static float laplacian(float *x, int w, int h, int i, int j)
{
	getpixel_operator p = getpixel_1;

	float r = -4 * p(x, w, h, i  , j  )
		     + p(x, w, h, i+1, j  )
		     + p(x, w, h, i  , j+1)
		     + p(x, w, h, i-1, j  )
		     + p(x, w, h, i  , j-1);

	return r;
}

// evaluate the laplacian of image x at point i, j
static float bilaplacian(float *x, int w, int h, int i, int j)
{
	float r = -4 * laplacian(x, w, h, i  , j  )
		     + laplacian(x, w, h, i+1, j  )
		     + laplacian(x, w, h, i  , j+1)
		     + laplacian(x, w, h, i-1, j  )
		     + laplacian(x, w, h, i  , j-1);

	return r;
}

// returns the largest change performed all over the image
static float perform_one_iteration(float *x, int w, int h,
		int (*mask)[2], int nmask, float tstep)
{
	float maxupdate = 0;
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i;

		float new = x[idx] + tstep * bilaplacian(x, w, h, i, j);

		float update = fabs(x[idx] - new);
		if (update > maxupdate)
			maxupdate = update;

		x[idx] = new;
	}
	return maxupdate;
}

// build a mask of the NAN positions on image "x"
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[2]
{
	int nmask = 0;
	for (int i = 0; i < w*h; i++)
		if (isnan(x[i]))
			nmask += 1;
	int (*mask)[2] = xmalloc(w*h*2*sizeof(int)), cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isnan(x[j*w + i])) {
			mask[cx][0] = i;
			mask[cx][1] = j;
			cx += 1;
		}
	assert(cx == nmask);

	*out_nmask = nmask;
	return mask;
}

// fill the holes of the image x using an harmonic function
static void harmonic_extension_with_init(
		float *y,        // output image
		float *x,        // input image (NAN values indicate holes)
		int w,           // image width
		int h,           // image height
		float timestep,  // time step for the numerical scheme
		int niter,       // number of iterations to run
		float *initialization
		)
{
	// build list of masked pixels
	int nmask, (*mask)[2] = build_mask(&nmask, x, w, h);

	// initialize the solution to the given data at the masked pixels
	for (int i = 0; i < w*h; i++)
		y[i] = isfinite(x[i]) ? x[i] : initialization[i];

	// do the requested iterations
	for (int i = 0; i < niter; i++)
	{
		float u = perform_one_iteration(y, w, h, mask, nmask, timestep);

		if (u < 1e-10) break;

		//if (0 == i % 10)
		//fprintf(stderr, "size = %dx%d, iter = %d, maxupdate = %g\n",
		//		w, h, i, u);
	}

	free(mask);
}

#include <stdbool.h>
#include "smapa.h"
SMART_PARAMETER(PREFILTER,0)

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	if (PREFILTER() > 0) {
		float *fin = xmalloc(iw*ih*sizeof*fin);
		for (int j = 0; j < ih; j++)
		for (int i = 0; i < iw; i++)
		{
			//if (isfinite(in[iw*j+i])) {
			//	float a[9], m = 0;
			//	for (int ii=-1;ii<=1;ii++)
			//	for (int jj=-1;jj<=1;jj++)
			//		a[3*ii+jj] = p(in, iw, ih, i+ii, j+jj);
			//	int cx = 0;
			//	for (int k = 0; k < 9; k++)
			//		if (isfinite(a[k])) {
			//			m += a[k];
			//			cx += 1;
			//		}
			//	fin[iw*j + i] = cx ? m/cx : NAN;
			//} else {
				float a[4], m = 0;
				a[0] = p(in, iw, ih, i+1, j);
				a[1] = p(in, iw, ih, i-1, j);
				a[2] = p(in, iw, ih, i, j+1);
				a[3] = p(in, iw, ih, i, j-1);
				int cx = 0;
				for (int k = 0; k < 4; k++)
					if (isfinite(a[k])) {
						m += a[k];
						cx += 1;
					}
				fin[iw*j + i] = cx ? m/cx : NAN;
			//}
		}
		in = fin;
	}
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, 2*i, 2*j);
		a[1] = p(in, iw, ih, 2*i+1, 2*j);
		a[2] = p(in, iw, ih, 2*i, 2*j+1);
		a[3] = p(in, iw, ih, 2*i+1, 2*j+1);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i] = cx ? m/cx : NAN;
	}
	if (PREFILTER() > 0) free(in);
}

// evaluate a bilinear cell at the given point
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// evaluate an image at a sub-pixel position, using bilinear interpolation
static float bilinear_interpolation(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpixel_1(x, w, h, ip  , iq  );
	float b = getpixel_1(x, w, h, ip+1, iq  );
	float c = getpixel_1(x, w, h, ip  , iq+1);
	float d = getpixel_1(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		out[ow*j+i] = p(in, iw, ih, round((i-0.5)/2), round((j-0.5)/2));
		//out[ow*j+i] = bilinear_interpolation(in, iw, ih, (i-0.5)/2, (j-0.5)/2);
}

// extension of an image by laplace equation
void elap_recursive(float *out, float *in, int w, int h,
		float timestep, int niter, int scale)
{
	float *init = xmalloc(w*h*sizeof*init);
	if (scale > 1)
	{
		int ws = ceil(w/2.0);
		int hs = ceil(h/2.0);
		float *ins  = xmalloc(ws * hs * sizeof*ins);
		float *outs = xmalloc(ws * hs * sizeof*outs);
		zoom_out_by_factor_two(ins, ws, hs, in, w, h);
		elap_recursive(outs, ins, ws, hs, timestep, niter, scale - 1);
		zoom_in_by_factor_two(init, w, h, outs, ws, hs);

		free(ins);
		free(outs);
	} else {
		for (int i = 0 ; i < w*h; i++)
			init[i] = 0;
	}
	harmonic_extension_with_init(out, in, w, h, timestep, niter, init);
	free(init);
}

// extension by laplace equation of each channel of a color image
void elap_recursive_separable(float *out, float *in, int w, int h, int pd,
		float timestep, int niter, int scale)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		elap_recursive(outl, inl, w, h, timestep, niter, scale);
	}
}

#define MAIN_ELAP_RECSEP

#ifdef MAIN_ELAP_RECSEP
#include "iio.h"

int main(int argc, char *argv[])
{
	if (argc != 7) {
		fprintf(stderr, "usage:\n\t"
		"%s TSTEP NITER NS data.png mask.png out.png\n", *argv);
		//0 1     2     3  4        5        6
		return 1;
	}
	float timestep = atof(argv[1]);
	int niter = atoi(argv[2]);
	int nscales = atoi(argv[3]);
	char *filename_in = argv[4];
	char *filename_mask = argv[5];
	char *filename_out = argv[6];

	int w[2], h[2], pd;
	float *in = iio_read_image_float_split(filename_in, w, h, &pd);
	float *mask = iio_read_image_float(filename_mask, w+1, h+1);
	if (w[0] != w[1] || h[0] != h[1])
		return fprintf(stderr, "image and mask file size mismatch"
				" %d %d != %d %d", w[0], h[0], w[1], h[1]);
	float *out = xmalloc(*w**h*pd*sizeof*out);

	for (int i = 0; i < *w * *h; i++)
		if (mask[i] > 0)
			for (int l = 0; l < pd; l++)
				in[*w**h*l+i] = NAN;

	elap_recursive_separable(out, in, *w, *h, pd, timestep, niter, nscales);

	iio_write_image_float_split(filename_out, out, *w, *h, pd);

	return 0;
}
#endif
