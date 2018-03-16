// tools for multi-scale optical flow {{{1
//
//
// There are two natural definitions for multi-scale optical flow:
//
// 1. Multi-scale by projection (can be used with any optical flow method)
// 2. Perturbative multi-scale (can be used with methods based on optimization)
//


// #includes {{{1
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "iio.h"

#include "smapa.h"

#include "fail.c"
#include "xmalloc.c"

// typedefs {{{1
typedef void (*generic_optical_flow)(
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		int width, int height,
		void *data);

typedef void (*iterative_optical_flow)(
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		float *in_u, float *in_v,
		int width, int height,
		void *data);

typedef float (*extension_operator_float)(
		float *image, int width, int height,
		int i, int j);

typedef float (*interpolation_operator_float)(
		float *image, int width, int height,
		float x, float y);

// utility functions {{{1
static float extend_float_image_constant(float *xx, int w, int h, int i, int j)
{
	float (*x)[w] = (void*)xx;
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j][i];
}

static float cell_interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float cell_interpolate_nearest(float a, float b, float c, float d,
					float x, float y)
{
	// return a;
	if (x<0.5) return y<0.5 ? a : b;
	else return y<0.5 ? c : d;
}

static float cell_interpolate(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	case 0: return cell_interpolate_nearest(a, b, c, d, x, y);
	//case 1: return marchi(a, b, c, d, x, y);
	case 2: return cell_interpolate_bilinear(a, b, c, d, x, y);
	default: return 0;
	}
}

static float interpolate_float_image_bilinearly(float *x, int w, int h,
		float i, float j)
{
	int ii = i;
	int jj = j;
	extension_operator_float p = extend_float_image_constant;
	float a = p(x, w, h, ii  , jj  );
	float b = p(x, w, h, ii  , jj+1);
	float c = p(x, w, h, ii+1, jj  );
	float d = p(x, w, h, ii+1, jj+1);
	return cell_interpolate(a, b, c, d, i-ii, j-jj, 2);
}

static void save_debug_image(char *fpat, int id, float *x, int w, int h)
{
	return;
	char filename[0x100];
	snprintf(filename, 0x100, fpat, id);
	fprintf(stderr, "saving image \"%s\"\n", filename);
	iio_write_image_float(filename, x, w, h);
}

static void save_debug_flow(char *fpat, int id, float*u, float*v, int w, int h)
{
	return;
	char filename[0x100];
	snprintf(filename, 0x100, fpat, id);
	float *f = xmalloc(w*h*2*sizeof*f);
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}

	fprintf(stderr, "saving field \"%s\"\n", filename);
	iio_write_image_float_vec(filename, f, w, h, 2);
	xfree(f);
}

// image scaling {{{1
//static void downscale_image_old(float *out, float *in,
//		int outw, int outh, int inw, int inh,
//		float scalestep)
//{
//	assert(scalestep == 2);
//	fprintf(stderr, "downscale(%g): %dx%d => %dx%d\n",
//			scalestep, inw, inh, outw, outh);
//	assert(2*outw >= inw);
//	assert(2*outh >= inh);
//
//	extension_operator_float p = extend_float_image_constant;
//
//	for (int j = 0; j < outh; j++)
//	for (int i = 0; i < outw; i++)
//		out[outw*j+i] = 0.25 * ( p(in,inw,inh, 2*i, 2*j)
//		       + p(in,inw,inh, 2*i+1, 2*j)
//		       + p(in,inw,inh, 2*i, 2*j+1)
//		       + p(in,inw,inh, 2*i+1, 2*j+1));
//}

static void downsa_v2(float *out, float *in,
		int outw, int outh, int inw, int inh)
{
	//fprintf(stderr, "\n\ndownsa_v2 (%d x %d) => (%d x %d)\n\n",
	//		inw,inh,outw,outh);
	if (2*outw != inw) fail("bad horizontal halving");
	if (2*outh != inh) fail("bad vertical halving");

	float (*y)[outw] = (void*)out;
	float (*x)[inw] = (void*)in;

	for (int j = 0; j < outh; j++)
	for (int i = 0; i < outw; i++) {
		float g = 0;
		g += x[2*j+0][2*i+0];
		g += x[2*j+0][2*i+1];
		g += x[2*j+1][2*i+0];
		g += x[2*j+1][2*i+1];
		y[j][i] = g/4;
	}
}

SMART_PARAMETER(MAGIC_SIGMA,1.6)
SMART_PARAMETER(PRESMOOTH,0)

static void downscale_image(float *out, float *in,
		int outw, int outh, int inw, int inh,
		float scalestep)
{
	if (scalestep == -2) {downsa_v2(out,in,outw,outh,inw,inh); return;}
	//fprintf(stderr, "downscale(%g): %dx%d => %dx%d\n",
	//		scalestep, inw, inh, outw, outh);

	assert(scalestep > 1);
	assert(scalestep * outw >= inw);
	//assert(scalestep * outw <= inw + 1);
	assert(scalestep * outh >= inh);
	//assert(scalestep * outh <= inh + 1);

	float factorx = inw/(float)outw;
	float factory = inh/(float)outh;

	float blur_size = MAGIC_SIGMA()*sqrt((factorx*factory-1)/3);

	fprintf(stderr, "blur_size = %g\n", blur_size);

	float *gin = xmalloc(inw * inh * sizeof(float));
	if (outw < inw || outh < inh) {
		void gblur_gray(float*, float*, int, int, float);
		gblur_gray(gin, in, inw, inh, blur_size);
	} else {
		assert(inw == outw);
		assert(inh == outh);
		for (int i = 0; i < inw*inh; i++)
			gin[i] = in[i];
	}

	// XXX ERROR FIXME
	// TODO: zoom by fourier, or zoom by bicubic interpolation
	interpolation_operator_float ev = interpolate_float_image_bilinearly;

	for (int j = 0; j < outh; j++)
	for (int i = 0; i < outw; i++)
	{
		float x = factorx*i;
		float y = factory*j;
		out[outw*j + i] = ev(gin, inw, inh, x, y);
	}

	xfree(gin);
}

// starting from a high-resolution image, produce a pyramid of lower-resolution
// versions
static void produce_upwards_pyramid(
		float **out_pyrx, int *out_pyrw, int *out_pyrh,
		float *x, int w, int h,
		int nscales, float sscalestep)
{
	float scalestep = fabs(sscalestep);
	assert(scalestep > 1);

	// compute pyramid sizes
	int pyrsize[nscales][2];
	pyrsize[0][0] = w;
	pyrsize[0][1] = h;
	for (int i = 1; i < nscales; i++) {
		pyrsize[i][0] = ceil(pyrsize[i-1][0] / scalestep);
		pyrsize[i][1] = ceil(pyrsize[i-1][1] / scalestep);
	}

	// save requested outputt
	for (int i = 0; i < nscales; i++) {
		if (out_pyrw) out_pyrw[i] = pyrsize[i][0];
		if (out_pyrh) out_pyrh[i] = pyrsize[i][1];
	}

	if (!out_pyrx) return;

	// alloc pyramid levels
	for (int i = 0; i < nscales; i++)
		out_pyrx[i] = xmalloc(pyrsize[i][0] * pyrsize[i][1] * sizeof*x);

	if (!x) return;

	// fill initial level
	//downscale_image(pyr[0], x, w, h, w, h, 1.0);
	for (int i = 0; i < w*h; i++)
		out_pyrx[0][i] = x[i];
	if (PRESMOOTH() > 0) {
		float presmooth = PRESMOOTH();
		void gblur_gray(float*, float*, int, int, float);
		gblur_gray(out_pyrx[0], x, w, h, presmooth);
	} else
		for (int i = 0; i < w*h; i++)
			out_pyrx[0][i] = x[i];


	// propagate information to other levels
	for (int i = 1; i < nscales; i++)
		downscale_image(out_pyrx[i], out_pyrx[i-1],
				pyrsize[i][0], pyrsize[i][1],
				pyrsize[i-1][0], pyrsize[i-1][1],
				sscalestep);

}

static void dealloc_pyramid(float **pyr, int n)
{
	for (int i = 0; i < n; i++)
		xfree(pyr[i]);
}

static void upscale_image(float *out, float *in,
		int outw, int outh, int inw, int inh,
		float sscalestep)
{
	float scalestep = fabs(sscalestep);
	assert(scalestep > 1);
	//assert(outw >= scalestep*inw);
	//assert(outw <= scalestep*(inw + 1));
	//assert(outh >= scalestep*inh);
	//assert(outh <= scalestep*(inh + 1));

	float factorx = outw/(float)inw;
	float factory = outh/(float)inh;

	interpolation_operator_float ev = interpolate_float_image_bilinearly;

	for (int j = 0; j < outh; j++)
	for (int i = 0; i < outw; i++)
	{
		float x = i/factorx;
		float y = j/factory;
		out[outw*j+i] = ev(in, inw, inh, x, y);
	}
}

static void upscale_field(
		float *outu, float *outv,
		float *inu, float *inv,
		int outw, int outh,
		int inw, int inh,
		float scalestep)
{
	upscale_image(outu, inu, outw, outh, inw, inh, scalestep);
	upscale_image(outv, inv, outw, outh, inw, inh, scalestep);

	float factorx = outw/(float)inw;
	float factory = outh/(float)inh;

	for (int i = 0; i < outw*outh; i++)
	{
		outu[i] *= factorx;
		outv[i] *= factory;
	}
}


//static void upscale_field_old(float *out, float *in,
//		int outw, int outh, int inw, int inh,
//		float scalestep)
//{
//	assert(scalestep == 2);
//	assert(outw <= 2*inw);
//	assert(outh <= 2*inh);
//
//	float (*o)[outw] = (void*)out;
//	float (*x)[inw] = (void*)in;
//
//	for (int j = 0; j < outh; j++)
//	for (int i = 0; i < outw; i++)
//		o[j][i] = 0;
//
//	// TODO: think about anchoring here, probably this is 0.5 pixels off!
//	for (int j = 0; j < inh; j++)
//	for (int i = 0; i < inw; i++)
//	{
//		o[2*j][2*i] = x[j][i];
//		if (2*j+1 < outh && j+1 < inh)
//			o[2*j+1][2*i] = (x[j][i] + x[j+1][i])/2;
//		if (2*i+1 < outw && i+1 < inw)
//			o[2*j][2*i+1] = (x[j][i] + x[j][i+1])/2;
//		if (2*j+1 < outh && 2*i+1 < outw && j+1 < inh && i+1 < inw)
//			o[2*j+1][2*i+1] = (x[j][i]+x[j][i+1]+x[j+1][i]+x[j+1][i+1])/4;
//	}
//
//	for (int i = 0; i < outw * outh; i++)
//		out[i] *= scalestep;
//}


// image warping {{{1
static void backwarp_image(float *x_out, float *x_in, float *u, float *v,
		int w, int h)
{
	interpolation_operator_float ev = interpolate_float_image_bilinearly;

	float *tmpx = xmalloc(w*h*sizeof*tmpx);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float p = i + u[w*j+i];
		float q = j + v[w*j+i];
		tmpx[w*j+i] = ev(x_in, w, h, p, q);
	}

	for (int i = 0; i < w*h; i++)
		x_out[i] = tmpx[i];

	xfree(tmpx);
}

// multiscale flow {{{1

// ungenericizer {{{2
static int global_idx;

static int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

SMART_PARAMETER(FLOW_MS_DO_FILTER,0)

static void filter_of(float *uu, float *vv, int w, int h)
{
	if (!FLOW_MS_DO_FILTER()) return;
	float *o_uu = xmalloc(w * h * sizeof*uu);
	float *o_vv = xmalloc(w * h * sizeof*uu);
	for (int i = 0; i < w*h; i++) {
		o_uu[i] = uu[i];
		o_vv[i] = vv[i];
	}
	float (*u)[w] = (void*)uu;
	float (*v)[w] = (void*)vv;
	float (*ou)[w] = (void*)o_uu;
	float (*ov)[w] = (void*)o_vv;
	for (int j = 1; j < h-1; j++)
	for (int i = 1; i < w-1; i++) {
		float mu[9], mv[9];
		int idx = 0;
		for (int di = -1; di <= 1; di++)
		for (int dj = -1; dj <= 1; dj++) {
			mu[idx] = u[j+dj][i+di];
			mv[idx] = v[j+dj][i+di];
			idx += 1;
		}
		qsort(mu, 9, sizeof*mu, compare_floats);
		qsort(mv, 9, sizeof*mu, compare_floats);
		ou[j][i] = mu[4];
		ov[j][i] = mv[4];
	}
	for (int i = 0; i < w*h; i++) {
		uu[i] = o_uu[i];
		vv[i] = o_vv[i];
	}
	free(o_uu);
	free(o_vv);
}

// this function should be a closure to turn a generic optical flow
// into an iterative optical flow
static void iteritized(generic_optical_flow of,
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		float *in_u, float *in_v,
		int w, int h,
		void *data)
{
	float *wb = xmalloc(w * h * sizeof*wb);
	float *u = xmalloc(w * h * sizeof*wb);
	float *v = xmalloc(w * h * sizeof*wb);
	float *wu = xmalloc(w * h * sizeof*wb);
	float *wv = xmalloc(w * h * sizeof*wb);

	backwarp_image(wb, in_b, in_u, in_v, w, h);

	save_debug_image("/tmp/ms_debug_image_wb0_%02d", global_idx, wb, w, h);

	of(u, v, in_a, wb, w, h, data);
	filter_of(u, v, w, h);
	for (int i = 0; i < w*h; i++)
	{
		float r = hypot(u[i], v[i]);
		float maxr = 1;
		if (r > maxr) {
			u[i] /= (r/maxr);
			v[i] /= (r/maxr);
		}
	}


	backwarp_image(wu, in_u, u, v, w, h);
	backwarp_image(wv, in_v, u, v, w, h);

	save_debug_flow("/tmp/ms_debug_field_duv_%02d", global_idx, u, v, w, h);
	for (int i = 0; i < w *h; i++) {
		out_u[i] = wu[i] + u[i];
		out_v[i] = wv[i] + v[i];
	}
	save_debug_flow("/tmp/ms_debug_field_uv_%02d", global_idx,
			out_u, out_v, w, h);
	xfree(wv);
	xfree(wu);
	xfree(v);
	xfree(u);
	xfree(wb);
}

// this function simply copies the flow
static void do_nothing(generic_optical_flow of,
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		float *in_u, float *in_v,
		int w, int h,
		void *data)
{
	for (int i = 0; i < w *h; i++) {
		out_u[i] = in_u[i];
		out_v[i] = in_v[i];
	}
	save_debug_flow("/tmp/ms_debug_field_uv_%02d", global_idx,
			out_u, out_v, w, h);
}

// perturbative multiscale {{{2
void perturbative_multi_scale_optical_flow(float *u, float *v,
		float *a, float *b, int w, int h,
		iterative_optical_flow of, void *data,
		int nscales, float scalestep)
{
	fail("not yet implemented");
}

SMART_PARAMETER_SILENT(NWARPS,1)

// generic multiscale {{{2
void generic_multi_scale_optical_flow(float *out_u, float *out_v,
		float *in_a, float *in_b, int in_w, int in_h,
		generic_optical_flow of, void *data,
		int nscales, float sstep, int last_scale)
{
	if (last_scale >= nscales) last_scale = nscales - 1;
	float step = fabs(sstep);
	assert(step > 1);

	float *a[nscales], *b[nscales], *u[nscales], *v[nscales];
	int w[nscales], h[nscales];

	//                      op ow oh ix
	produce_upwards_pyramid(0, w, h, 0,    in_w, in_h, nscales,sstep);
	produce_upwards_pyramid(a, 0, 0, in_a, in_w, in_h, nscales,sstep);
	produce_upwards_pyramid(b, 0, 0, in_b, in_w, in_h, nscales,sstep);
	produce_upwards_pyramid(u, 0, 0, 0,    in_w, in_h, nscales,sstep);
	produce_upwards_pyramid(v, 0, 0, 0,    in_w, in_h, nscales,sstep);

	for (int i = 0; i < nscales; i++) {
		save_debug_image("/tmp/ms_debug_image_a_%02d", i,
				a[i], w[i], h[i]);
		save_debug_image("/tmp/ms_debug_image_b_%02d", i,
				b[i], w[i], h[i]);
	}

	int nwarps = NWARPS();

	int s = nscales - 1;
	for (int i = 0; i < w[s]*h[s]; i++)
		u[s][i] = v[s][i] = 0;
	while (s >= 0) {
		global_idx = s;

		// run flow at this level
		if (s >= last_scale)
			for (int i = 0; i < nwarps; i++)
				iteritized(of, u[s], v[s], a[s], b[s],
						u[s], v[s], w[s], h[s], data);
		else
			do_nothing(of, u[s], v[s], a[s], b[s], u[s], v[s],
							w[s], h[s], data);

		if (!s) break;

		// upscale flow to the next level
		upscale_field(u[s-1], v[s-1], u[s], v[s],
		              w[s-1], h[s-1], w[s], h[s], step);

		save_debug_flow("/tmp/ms_debug_field_suv_%02d", s,
				u[s-1], v[s-1], w[s-1], h[s-1]);

		// iterate
		s = s - 1;
	}

	for (int i = 0; i < in_w * in_h; i++) {
		out_u[i] = u[0][i];
		out_v[i] = v[0][i];
	}

	dealloc_pyramid(a, nscales);
	dealloc_pyramid(b, nscales);
	dealloc_pyramid(u, nscales);
	dealloc_pyramid(v, nscales);
}

// genericized flows {{{1
static void genericized_hs(
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		int width, int height,
		void *data)
{
	float *fdata = data;
	float alpha = fdata[0];
	int niter = fdata[1];
	float epsilon = fdata[2];
	void hs(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha);
	void hs_stopping(float *u, float *v, float *a, float *b, int w, int h,
		int niter, float alpha, float epsilon);
	//fprintf(stderr, "calling HS with input of size %dx%d\n", width, height);
	//hs(out_u, out_v, in_a, in_b, width, height, niter, alpha);
	hs_stopping(out_u, out_v, in_a, in_b, width, height,
			niter, alpha, epsilon);
}

static void genericized_lk(
		float *out_u, float *out_v,
		float *in_a, float *in_b,
		int w, int h,
		void *data)
{
	float *fdata = data;
	float kside = fdata[0];
	float sigma = fdata[1];
	void least_squares_ofc(float *u, float *v,
		float *a, float *b, int w, int h, int kside, float sigma);
	least_squares_ofc(out_u, out_v, in_a, in_b, w, h, kside, sigma);
}

// actually usable api {{{1
void multi_scale_optical_flow(char *algorithm_name, float *pars, int npars,
		float *u, float *v, float *a, float *b, int w, int h,
		int nscales, float scale_step, int last_scale)
{
	fprintf(stderr, "multi-scale flow \"%s\" [%d %g %d]", algorithm_name,
			nscales, scale_step, last_scale);
	for (int i = 0; i < npars; i++)
		fprintf(stderr, " %g", pars[i]);
	fprintf(stderr, "\n");

	generic_optical_flow of = NULL;
	if (0 == strcmp("hs", algorithm_name)) {
		if (npars != 3)
			fail("flow \"%s\" needs 3 parameters", algorithm_name);
		of = genericized_hs;
	}
	if (0 == strcmp("lk", algorithm_name)) {
		if (npars != 2)
			fail("flow \"%s\" needs 2 parameters", algorithm_name);
		of = genericized_lk;
	}
	if (0 == strcmp("zero", algorithm_name)) {
		for (int i = 0; i < w*h; i++)
			u[i] = v[i] = 0;
		return;
	}
	if (!of) fail("unrecognized flow method \"%s\"", algorithm_name);

	generic_multi_scale_optical_flow(u, v, a, b, w, h, of, pars,
			nscales, scale_step, last_scale);
}

void multi_scale_optical_flow_pd(char *algorithm_name, float *pars, int npars,
		float *u, float *v, float *a, float *b, int w, int h, int pd,
		int nscales, float scale_step, int last_scale)
{
	float *ga = xmalloc(w*h*sizeof*ga);
	float *gb = xmalloc(w*h*sizeof*ga);
	for (int i = 0; i < w*h; i++)
	{
		float ma = 0;
		float mb = 0;
		for (int j = 0; j < pd; j++) {
			ma += a[i*pd + j];
			mb += b[i*pd + j];
		}
		ga[i] = ma/pd;
		gb[i] = mb/pd;
	}

	multi_scale_optical_flow(algorithm_name, pars, npars,
			u, v, ga, gb, w, h, nscales, scale_step, last_scale);

	xfree(ga);
	xfree(gb);
}

// main for testing api {{{1
#ifdef USE_MAINAPI

#define BAD_MIN(a,b) (a)<(b)?(a):(b)

static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

int main(int argc, char *argv[])
{
	if (argc != 9 && argc != 8) {
		fprintf(stderr, "usage:\n\t"
		"%s a b method \"params\" step nscales lastscale [f]\n", *argv);
	//       0  1 2 3        4        5    6       7          8
		return EXIT_FAILURE;
	}
	char *filename_a = argv[1];
	char *filename_b = argv[2];
	char *method_id = argv[3];
	char *parstring = argv[4];
	float scale_step = atof(argv[5]);
	int nscales = atoi(argv[6]);
	int last_scale = atoi(argv[7]);
	char *filename_f = argc > 8 ? argv[8] : "-";

	float params[0x100];
	int nparams = parse_floats(params, 0x100, parstring);

	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh) fail("input image size mismatch");
	float *u = xmalloc(2 * w * h * sizeof(float));
	float *v = u + w * h;

	// bound the number of scales to disallow images smaller than 16 pixels
	float Nscales = 1.5+log(BAD_MIN(w,h)/3.0)/log(fabs(scale_step));
	if (Nscales < nscales)
		nscales = Nscales;

	multi_scale_optical_flow(method_id, params, nparams,
			u, v, a, b, w, h, nscales, scale_step, last_scale);

	iio_write_image_float_split(filename_f, u, w, h, 2);

	xfree(u);
	xfree(a);
	xfree(b);
	return EXIT_SUCCESS;
}
#endif//USE_MAINAPI

// main for testing pyramid construction {{{1
#ifdef USE_MAINPCS

#include "iio.h"

int main(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr, "usage:\n\t%s img step nscales outpat\n", *v);
		//                          0 1   2    3       4
		return EXIT_FAILURE;
	}
	char *filename_in = v[1];
	float scalestep = atof(v[2]);
	int nscales = atoi(v[3]);
	char *filepattern_out = v[4];

	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	float *pyrx[nscales];
	int pyrw[nscales], pyrh[nscales];

	produce_upwards_pyramid(pyrx, pyrw, pyrh,
			x, w, h, nscales, scalestep);

	for (int i = 0; i < nscales; i++) {
		char buf[0x100];
		snprintf(buf, 0x100, filepattern_out, i);
		printf("SCALE NUMBER %d: %dx%d (%d) => \"%s\"\n", i,
				pyrw[i], pyrh[i], buf);
		iio_write_image_float(buf, pyrx[i], pyrw[i], pyrh[i]);
	}


	return EXIT_SUCCESS;
}

#endif//USE_MAIN


// main for testing generic multiscale {{{1
#ifdef USE_MAINPHS


#define BAD_MIN(a,b) (a)<(b)?(a):(b)

int main(int argc, char *argv[])
{
	if (argc != 9) {
		fprintf(stderr, "usage:\n\t"
			"%s a b alpha niter eps step nscales f\n", *argv);
		//       0  1 2 3     4     5   6    7       8
		return EXIT_FAILURE;
	}
	char *filename_a = argv[1];
	char *filename_b = argv[2];
	float alpha = atof(argv[3]);
	int niter = atoi(argv[4]);
	float epsilon = atof(argv[5]);
	float scalestep = atof(argv[6]);
	int nscales = atoi(argv[7]);
	char *filename_f = argv[8];

	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh) fail("input image size mismatch");
	float *u = xmalloc(w * h * sizeof(float));
	float *v = xmalloc(w * h * sizeof(float));

	float Nscales = 1.5+log(BAD_MIN(w,h)/3.0)/log(fabs(scalestep));
	if (Nscales < nscales)
		nscales = Nscales;

	float fdata[3] = {alpha, niter, epsilon};
	generic_multi_scale_optical_flow(u, v, a, b, w, h,
			genericized_lk, fdata, nscales, scalestep, 0);

	float *f = xmalloc(w*h*2*sizeof*f);
	for (int i = 0; i < w*h; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_write_image_float_vec(filename_f, f, w, h, 2);

	xfree(f);
	xfree(u);
	xfree(v);
	xfree(a);
	xfree(b);

	return EXIT_SUCCESS;
}
#endif//USE_MAINPHS

// vim:set foldmethod=marker:
