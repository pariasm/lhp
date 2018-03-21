// various program to filter segments

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

int main_length_above_threshold(int c, char *v[])
{
	if (c != 2) return fprintf(stderr, "usage:\n\t"
			"%s dist < in.txt > out.txt\n", *v);
	//                0 1
	float t = atof(v[1]);

	float x[4];
	while (4 == scanf("%g %g %g %g\n", x+0, x+1, x+2, x+3))
		if (hypot(x[0] - x[2], x[1] - x[3]) > t)
			printf("%g %g %g %g\n", x[0], x[1], x[2], x[3]);

	return 0;
}

static void apply_homography(float *y, double H[9], float *x)
{
	double p = H[0] * x[0] + H[1] * x[1] + H[2];
	double q = H[3] * x[0] + H[4] * x[1] + H[5];
	double r = H[6] * x[0] + H[7] * x[1] + H[8];
	y[0] = p / r;
	y[1] = q / r;
}

#include "parsenumbers.c"
int main_transform_by_homography(int c, char *v[])
{
	if (c != 2) return fprintf(stderr, "usage:\n\t"
			"%s homography < in.txt > out.txt\n", *v);
	//                0 1
	double H[9];
	read_n_doubles_from_string(H, v[1], 9);

	float x[4];
	while (4 == scanf("%g %g %g %g\n", x+0, x+1, x+2, x+3))
	{
		float y[4];
		apply_homography(y+0, H, x+0);
		apply_homography(y+2, H, x+2);
		printf("%g %g %g %g\n", y[0], y[1], y[2], y[3]);
	}

	return 0;
}


#include "iio.h"
int main_completely_inside_binary_mask(int c, char *v[])
{
	if (c != 2) return fprintf(stderr, "usage:\n\t"
			"%s mask.png < in.txt > out.txt\n", *v);
	//                0 1
	int w, h;
	uint8_t *m = iio_read_image_uint8(v[1], &w, &h);


	float x[4];
	while (4 == scanf("%g %g %g %g\n", x+0, x+1, x+2, x+3))
	//                  0  1  2  3
	{
		int p[4] = { x[0], x[1], x[2], x[3] };
		if (p[0] < 0 || p[0] >= w) continue;
		if (p[2] < 0 || p[2] >= w) continue;
		if (p[1] < 0 || p[1] >= h) continue;
		if (p[3] < 0 || p[3] >= h) continue;
		if (m[w*p[1] + p[0]] || m[w*p[3] + p[0]])
			printf("%g %g %g %g\n", x[0], x[1], x[2], x[3]);
	}

	return 0;
}

int main(int c, char *v[])
{
	if (c >= 2 && 0 == strcmp(v[1], "minlength"))
		return main_length_above_threshold(c-1, v+1);
	if (c >= 2 && 0 == strcmp(v[1], "inmask"))
		return main_completely_inside_binary_mask(c-1, v+1);
	if (c >= 2 && 0 == strcmp(v[1], "homography"))
		return main_transform_by_homography(c-1, v+1);
	return fprintf(stderr, "usage:\n\t"
			"segfilter {minlength|inmask|homography} ...\n");
}
