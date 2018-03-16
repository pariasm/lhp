// FLOWH: compute and visualize vector field histograms

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"

static int closed_bound(int a, int x, int b)
{
	if (x < a) return a;
	if (x <= b) return x;
	return b;
}

static bool unmasked(int idx)
{
	static bool initialized = false;
	static int w, h;
	static float *x = NULL;
	if (!initialized)
	{
		initialized = true;
		char *fname = getenv("FLOWH_MASK");
		if (!fname) goto endif;
		x = iio_read_image_float(fname, &w, &h);
		if (x) fprintf(stderr, "using mask file \"%s\"\n", fname);
	}
endif:
	if (x) {
		assert(idx >= 0);
		assert(idx < w*h);
	}
	return !x || (x[idx] > 0);
}

void flowh(float *y, int ow, int oh, float *x, int w, int h, int nbu, int nu)
{
	assert(ow == nbu * (2*nu + 1));
	assert(oh == ow);

	for (int i = 0; i < ow*oh; i++)
		y[i] = 0;

	int offs = ow/2;
	for (int i = 0; i < w*h; i++)
	{
		float *f = x + 2*i;
		int fi[2];
		for (int j = 0; j < 2; j++) {
			fi[j] = closed_bound(0, offs + round(f[j]*nbu), ow-1);
		}
		int idx = fi[1]*ow + fi[0];
		assert(idx >= 0);
		assert(idx < ow*oh);
		if (unmasked(i))
			y[idx] += 1;
	}
}

#ifndef OMIT_MAIN
#include "pickopt.c"
int main(int c, char *v[])
{
	bool use_logscale = pick_option(&c, &v, "l", NULL);
	if (c != 3 && c != 4 && c != 5) {
		fprintf(stderr, "usage:\n\t%s nbu nu [in [out]]\n", *v);
		//                          0 1   2   3   4
		return 1;
	}
	int number_of_bins_per_unit = atoi(v[1]);
	int number_of_units = atoi(v[2]);
	char *infile = c > 3 ? v[3] : "-";
	char *outfile = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(infile, &w, &h, &pd);
	if (pd != 2) fail("2D vector field expected");

	int ow = number_of_bins_per_unit * (2 * number_of_units + 1);
	int oh = ow;
	float *y = xmalloc(ow*oh*sizeof*y);
	flowh(y, ow, oh, x, w, h, number_of_bins_per_unit, number_of_units);
	if (use_logscale) for (int i = 0; i < ow*oh; i++) y[i] = log(1+y[i]);
	iio_write_image_float_vec(outfile, y, ow, oh, 1);
	return 0;
}
#endif
