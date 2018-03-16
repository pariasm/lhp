
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>



#include "iio.h"

#include "fail.c"
#include "xfopen.c"
#include "xmalloc.c"

#include "smapa.h"
SMART_PARAMETER(OVERCROSS,1)


int main(int c, char *v[])
{
	int w, h;
	uint8_t *img_in = iio_read_image_uint8(v[1], &w, &h);
	uint8_t *img_out_raw = xmalloc(3*w*h);
	uint8_t (*img_out)[w][3] = (void*)img_out_raw;
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
			for (int l = 0; l < 3; l++)
				img_out[j][i][l] = img_in[w*j+i];
	float p[2];
	while (2 == fscanf(stdin, "%g %g\n", p, p+1)) {
		int iox = p[0];
		int ioy = p[1];
		if (iox <= 0) iox = 1;
		if (ioy <= 0) ioy = 1;
		if (iox >= w-1) iox = w-2;
		if (ioy >= h-1) ioy = h-2;
		img_out[ioy][iox][0] = 255;
		img_out[ioy][iox][1] = 0;
		img_out[ioy][iox][2] = 0;
		if (!OVERCROSS()) continue;
		img_out[ioy-1][iox][0] = 255;
		img_out[ioy-1][iox][1] = 0;
		img_out[ioy-1][iox][2] = 0;
		img_out[ioy+1][iox][0] = 255;
		img_out[ioy+1][iox][1] = 0;
		img_out[ioy+1][iox][2] = 0;
		img_out[ioy][iox+1][0] = 255;
		img_out[ioy][iox+1][1] = 0;
		img_out[ioy][iox+1][2] = 0;
		img_out[ioy][iox-1][0] = 255;
		img_out[ioy][iox-1][1] = 0;
		img_out[ioy][iox-1][2] = 0;
	}
	fprintf(stderr, "ok, saving file\n");
	iio_write_image_uint8_vec("-", img_out_raw, w, h, 3);
	return 0;
}
