#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "getpixel.c"

static int bound(int x, int min, int max)
{
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

static void croparound(float *out, int cw, int ch,
		float *in, int w, int h, int pd, int x, int y)
{
	int x0 = bound(x - cw/2, 0, w - cw - 2);
	int y0 = bound(y - ch/2, 0, h - ch - 2);
	int xf = x0 + cw + 1;  if (xf > w-1) {xf = w-1;}
	int yf = y0 + ch + 1;  if (yf > h-1) {yf = h-1;}

	if (w < cw) {x0 = 0; xf = cw-1;}
	if (h < ch) {y0 = 0; yf = ch-1;}

	//assert(cw == xf - x0 - 1);
	//assert(ch == yf - y0 - 1);

	for (int j = 0; j < ch; j++)
	for (int i = 0; i < cw; i++)
	for (int l = 0; l < pd; l++)
	{
		int ii = x0 + i;
		int jj = y0 + j;
		float g = getsample_0(in, w, h, pd, ii, jj, l);
		setsample_0(out, cw, ch, pd, i, j, l, g);
	}
}


#include "iio.h"
#include "xmalloc.c"

// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
static char *pick_option(int *c, char ***v, char *o, char *d)
{
	int argc = *c;
	char **argv = *v;
	for (int i = 0; i < argc - 1; i++)
		if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
		{
			char *r = argv[i+1];
			*c -= 2;
			for (int j = i; j < argc - 1; j++)
				(*v)[j] = (*v)[j+2];
			return r;
		}
	return d;
}

int main(int c, char *v[])
{
	int relative = atoi(pick_option(&c, &v, "r", "0"));
	if (c < 5 || c > 7) {
		fprintf(stderr, "usage:\n\t%s x y w h [in [out]]\n", *v);
		//                          0 1 2 3 4  5   6
		return EXIT_FAILURE;
	}
	float x = atof(v[1]);
	float y = atof(v[2]);
	int cw = atoi(v[3]);
	int ch = atoi(v[4]);
	char *filename_in = c > 5 ? v[5] : "-";
	char *filename_out = c > 6 ? v[6] : "-";

	int w, h, pd;
	float *image_in = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *image_out = xmalloc(cw*ch*pd*sizeof*image_out);

	if (relative) {
		x *= w;
		y *= h;
	}
	croparound(image_out, cw, ch, image_in, w, h, pd, x, y);

	iio_write_image_float_vec(filename_out, image_out, cw, ch, pd);
	return EXIT_SUCCESS;
}
