// the simplest implementation of morphological operators


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "iio.h"


// utility function that returns a valid pointer to memory
static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new) {
		fprintf(stderr, "ERROR: out of memory when requesting "
			       "%zu bytes\n", size);
		abort();
	}
	return new;
}

// the type of the "getpixel" operators
typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by nan
static float getpixel_nan(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return NAN;
	return x[i + j*w];
}

// the type of a morphological operator
typedef	void (*morphological_operator)(float*,float*,int,int,int*);

// a structuring element is a list E of integers
// E[0] = number of pixels
// E[1] = 0 (flags, not used yet)
// (E[2], E[3]) = position of the center
// (E[4], E[5]) = first pixel
// (E[6], E[7]) = second pixel
// ...

void morsi_erosion(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fmin(a, p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

void morsi_dilation(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = -INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fmax(a, p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

static float median(float *a, int n)
{
	if (n < 1) return NAN;
	if (n == 1) return *a;
	if (n == 2) return (a[0] + a[1])/2;
	qsort(a, n, sizeof*a, compare_floats);
	if (0 == n%2)
		return (a[n/2]+a[1+n/2])/2;
	else
		return a[n/2];
}

void morsi_median(float *y, float *x, int w, int h, int *e)
{
	getpixel_operator p = getpixel_nan;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a[e[0]];
		int cx = 0;
		for (int k = 0; k < e[0]; k++)
		{
			float v = p(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]);
			if (isfinite(v))
				a[cx++] = v;
		}
		y[j*w+i] = median(a, cx);
	}
}

void morsi_opening(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	morsi_dilation(y, t, w, h, e);
	free(t);
}

void morsi_closing(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_dilation(t, x, w, h, e);
	morsi_erosion(y, t, w, h, e);
	free(t);
}

void morsi_gradient(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_erosion(a, x, w, h, e);
	morsi_dilation(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = b[i] - a[i];
	free(a);
	free(b);
}

void morsi_igradient(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_erosion(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] - t[i];
	free(t);
}

void morsi_egradient(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_dilation(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = t[i] - x[i];
	free(t);
}

void morsi_laplacian(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_erosion(a, x, w, h, e);
	morsi_dilation(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = (a[i] + b[i] - 2*x[i])/2;
	free(a);
	free(b);
}

void morsi_enhance(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_laplacian(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] - t[i];
	free(t);
}

void morsi_strange(float *y, float *x, int w, int h, int *e)
{
	float *a = xmalloc(w*h*sizeof*a);
	float *b = xmalloc(w*h*sizeof*b);
	morsi_opening(a, x, w, h, e);
	morsi_closing(b, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = b[i] - a[i];
	free(a);
	free(b);
}

void morsi_tophat(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_opening(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = x[i] - t[i];
	free(t);
}

void morsi_bothat(float *y, float *x, int w, int h, int *e)
{
	float *t = xmalloc(w*h*sizeof*t);
	morsi_closing(t, x, w, h, e);
	for (int i = 0; i < w*h; i++)
		y[i] = t[i] - x[i];
	free(t);
}

static int *build_disk(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4, elen = 2*side*side+4;
	int *e = xmalloc(elen*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
	for (int j = -radius-1; j <= radius+1; j++)
		if (hypot(i,j) < radius) {
			e[2*cx+4] = i;
			e[2*cx+5] = j;
			cx += 1;
		}
	assert(cx < side*side);
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

static int *build_hrec(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4;
	int *e = xmalloc((2*side+4)*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
		if (abs(i) < radius) {
			e[2*cx+4] = i;
			e[2*cx+5] = 0;
			cx += 1;
		}
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

static int *build_vrec(float radius)
{
	if (!(radius >1)) return NULL;
	int side = 2*radius+4;
	int *e = xmalloc((2*side+4)*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
		if (abs(i) < radius) {
			e[2*cx+4] = 0;
			e[2*cx+5] = i;
			cx += 1;
		}
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}


static int *build_structuring_element_from_string(char *name)
{
	static int cross[] = {5,0,  0,0, -1,0, 0,0, 1,0, 0,-1, 0,1 };
	static int square[] = {9,0, 0,0,
			-1,-1,-1,0,-1,1, 0,-1,0,0,0,1, 1,-1,1,0,1,1};
	int *r = NULL;
	if (0 == strcmp(name, "cross" )) r = cross;
	if (0 == strcmp(name, "square")) r = square;
	if (4 == strspn(name, "disk"  )) r = build_disk(atof(name+4));
	if (4 == strspn(name, "hrec"  )) r = build_hrec(atof(name+4));
	if (4 == strspn(name, "vrec"  )) r = build_vrec(atof(name+4));
	return r;
}


#include "iio.h"

static void demo(char *filename_out,
		morphological_operator op,
		int *structuring_element,
		float *x, int w, int h, int pd)
{
	// keep a counter to sort the images
	static int cx = 1;
	char buf[FILENAME_MAX];
	snprintf(buf, FILENAME_MAX, "%03d_%s", cx++, filename_out);

	// allocate space for the image
	float *y = malloc(w*h*pd*sizeof*y);

	// compute
	if (pd == 1)
		op(y, x, w, h, structuring_element);
	else
		for (int k = 0; k < pd; k++)
			op(y+k*w*h, x+k*w*h, w, h, structuring_element);

	// normalize the image, only in the laplacian case
	if (op == morsi_laplacian)
		for (int i = 0; i < w*h*pd; i++)
			y[i] = y[i] + 127;

	// save result
	iio_write_image_float_split(buf, y, w, h, pd);

	// cleanup
	free(y);
}

int main(int c, char **v)
{
	// process input arguments
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s element image\n", *v);
		//                          0 1       2
	int *structuring_element = build_structuring_element_from_string(v[1]);
	if (!structuring_element)
		return fprintf(stderr, "elements = cross, square, disk5 ...\n");
	char *filename_in = v[2];

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	// run the demo for each operation
	demo("erosion.png",   morsi_erosion,   structuring_element, x, w,h,pd);
	demo("dilation.png",  morsi_dilation,  structuring_element, x, w,h,pd);
	demo("median.png",    morsi_median,    structuring_element, x, w,h,pd);
	demo("opening.png",   morsi_opening,   structuring_element, x, w,h,pd);
	demo("closing.png",   morsi_closing,   structuring_element, x, w,h,pd);
	demo("gradient.png",  morsi_gradient,  structuring_element, x, w,h,pd);
	demo("igradient.png", morsi_igradient, structuring_element, x, w,h,pd);
	demo("egradient.png", morsi_egradient, structuring_element, x, w,h,pd);
	demo("laplacian.png", morsi_laplacian, structuring_element, x, w,h,pd);
	demo("enhance.png",   morsi_enhance,   structuring_element, x, w,h,pd);
	demo("tophat.png",    morsi_tophat,    structuring_element, x, w,h,pd);
	demo("bothat.png",    morsi_bothat,    structuring_element, x, w,h,pd);

	// cleanup and exit
	free(x);
	return 0;
}
