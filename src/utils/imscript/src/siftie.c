#include <assert.h>
#include <stdio.h>
#include <math.h>


#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"

#define SIFT_LENGTH 128
struct sift_keypoint {
	float pos[2], scale, orientation;
	float affinity[6];
	float id;
	float sift[SIFT_LENGTH];
};

struct ann_pair {
	int from, to;
	float v[2];
};

struct ann_trip {
	int froma, tob, toc;
	float v[3];
};

#include "smapa.h"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)



#if USE_KDTREE
#include <kdtree.h>
#endif

static void write_raw_sift(FILE *f, struct sift_keypoint *k)
{
	fprintf(f, "%g %g %g %g", k->pos[0], k->pos[1], k->scale,
			k->orientation);
	FORJ(SIFT_LENGTH)
		fprintf(f, " %g", k->sift[j]);
	fprintf(f, "\n");
}

//void write_raw_sift_bin(FILE *f, struct sift_keypoint *k)
//{
//	assert(SIFT_LENGTH==128);
//	size_t n = 4*sizeof(float)+128;
//	char buf[n];
//	float *p = (float*)buf;
//	p[0] = k->pos[0];
//	p[1] = k->pos[1];
//	p[2] = k->scale;
//	p[3] = k->orientation;
//	char *pp = (char *)(4+p);
//	FORI(128) pp[i] = k->sift[i];
//	size_t r = fwrite(buf, n, 1, f);
//	if (r != n) fail("could not write some sift descriptor");
//}
//
//SMART_PARAMETER_SILENT(SIFT_OUTBIN,0)
//SMART_PARAMETER_SILENT(SIFT_INBIN,0)

static void write_raw_sifts(FILE *f, struct sift_keypoint *k, int n)
{
	FORI(n) {
//		if (BIN_OSIFT() > 0)
//			write_raw_sift_bin(f, k);
//		else
			write_raw_sift(f, k);
		k += 1;
	}
}

static struct sift_keypoint *read_raw_sifts(FILE *f, int *no)
{
	int n, N = 132;
	float *t = read_ascii_floats(f, &n);
	if (n == 0) {*no=0; return NULL;}
	if (0 != n % N) fail("bad raw SIFT keypoints format (%d %d)", n, n%N);
	n /= N;
	struct sift_keypoint *r = xmalloc(n * sizeof * r);
	FORI(n) {
		int off = i*N;
		r->pos[0] = t[0+off];
		r->pos[1] = t[1+off];
		r->scale = t[2+off];
		r->orientation = t[3+off];
		FORJ(SIFT_LENGTH)
		{
			int c = t[4+j+off];
			if (c < 0 || c > 255) fail("bad sift value %d\n", c);
			r->sift[j] = c;
		}
		r += 1;
	}
	xfree(t);
	*no = n;
	return r-n;
}

static struct sift_keypoint *read_raw_sifts_fname(char *fname, int *n)
{
	FILE *f = xfopen(fname, "r");
	struct sift_keypoint *r = read_raw_sifts(f, n);
	xfclose(f);
	return r;
}

static void write_raw_siftb(FILE *f, struct sift_keypoint *k)
{
	float t[4] = {k->pos[0], k->pos[1], k->scale, k->orientation};
	size_t r = fwrite(t, sizeof(*t), 4, f);
	if (r != sizeof*t)
		fail("could not write SIFT descriptor table (r = %zu)", r);
	char T[SIFT_LENGTH];
	FORI(SIFT_LENGTH)
		if (k->sift[i] < 0 || k->sift[i] > 255)
			fail("bad sift descriptor entry %g", k->sift[i]);
		else
			T[i] = k->sift[i];
	r = fwrite(T, SIFT_LENGTH, 1, f);
	if (r != 1)
		fail("could not write SIFT descriptor table (r=%zu)", r);
}

static void write_raw_siftsb(FILE *f, struct sift_keypoint *k, int n)
{
	FORI(n)
		write_raw_siftb(f, k+i);
}

static void *freadwhole_f(FILE *f, long *on)
{
#if 0
	int r = fseek(f, 0, SEEK_END);
	if (r)
		error("can not determine size of file (%d)", r);
	long n = ftell(f);
	if (n < 0)
		error("can not determine size of file (%ld)", n);
	void *ret = xmalloc(n);
	long rr = fread(ret, 1, n, f);
	if (rr != n)
		error("could not read %ld bytes from file (only %ld)", n, rr);
	*on = n;
	return ret;
#else
	int r, n = 0, nt = 0;
	char *t = NULL;
	while(1) {
		if (n >= nt)
		{
			nt = 1000+2 * (nt + 1);
			t = xrealloc(t, nt);
		}
		r = fgetc(f);
		if (r == EOF)
			break;
		t[n] = r;
		n += 1;
	}
	*on = n;
	return t;
#endif
}

static void *freadwhole(const char *fname, long *on)
{
	FILE *f = xfopen(fname, "r");
	void *ret = freadwhole_f(f, on);
	xfclose(f);
	return ret;
}


// sift descriptor saved length
static const int SDSLEN = SIFT_LENGTH*1 + 4*sizeof(float);

static struct sift_keypoint *read_raw_siftsb(FILE *f, int *no)
{
	long nn;
	void *p = freadwhole_f(f, &nn);
	long n = nn / SDSLEN;
	if (n*SDSLEN != nn)
		fail("can not read binary sift file (%ld %ld)", n*SDSLEN, nn);
	struct sift_keypoint *ret = xmalloc(n * sizeof * ret);
	FORI(n) {
		struct sift_keypoint *k = ret + i;
		void *pi = SDSLEN * i + (char *)p;
		float *t = pi;
		k->pos[0] = t[0];
		k->pos[1] = t[1];
		k->scale = t[2];
		k->orientation = t[3];
		unsigned char *tt = 4*sizeof(float) + (unsigned char *)pi;
		FORJ(SIFT_LENGTH)
			k->sift[j] = tt[j];

	}
	xfree(p);
	*no = n;
	return ret;
}

SMART_PARAMETER(SIFT_BINARY,0)

static struct sift_keypoint *read_raw_sifts_gen(FILE *f, int *no)
{
	return SIFT_BINARY() ? read_raw_siftsb(f, no) : read_raw_sifts(f, no);
}
static void write_raw_sifts_gen(FILE *f, struct sift_keypoint *k, int n)
{
	SIFT_BINARY() ? write_raw_siftsb(f, k, n) : write_raw_sifts(f, k, n);
}
static void write_raw_sift_gen(FILE *f, struct sift_keypoint *k)
{
	SIFT_BINARY() ? write_raw_siftb(f, k) : write_raw_sift(f, k);
}

static float emvdistppf(float *a, float *b, int n)
{
	//float ac = a[0];
	//float bc = b[0];
	//float r = fabs(ac - bc);
	//FORI(n-1)
	//	r += fabs(ac + a[i+1] - bc - b[i+1]);
	//return r;
	float ac[n]; ac[0] = a[0]; FORI(n-1) ac[i+1] = ac[i] + a[i+1];
	float bc[n]; bc[0] = b[0]; FORI(n-1) bc[i+1] = bc[i] + b[i+1];
	float r = 0;
	FORI(n) r += fabs(ac[i] - bc[i]);
	fprintf(stderr, "%g\n", r);
	return r;
}

static double sqr(double x)
{
	return x*x;
}

static float distppft(float *a, float *b, int n, float t)
{
	double tt = t * t;
	double r = 0;
	for (int i = 0; i < n; i++)
		if (r > tt)
			return t + 1;
		else
			r = hypot(r, b[i] - a[i]);
			//r += sqr(b[i] - a[i]);
	return r;//sqrt(r);
}

static float distppf(float *a, float *b, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += sqr(b[i] - a[i]);
	return sqrt(r);
}

static float distlpf(float *a, float *b, int n, float p)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += pow(fabs(b[i] - a[i]), p);
	return pow(r, 1.0/p);
}

SMART_PARAMETER(DIST_DESCS_EMV,0)
SMART_PARAMETER(DIST_DESCS_LP,2)

static double dist_descs(struct sift_keypoint *a, struct sift_keypoint *b)
{
	//float af[SIFT_LENGTH]; FORI(SIFT_LENGTH) af[i] = a->sift[i];
	//float bf[SIFT_LENGTH]; FORI(SIFT_LENGTH) bf[i] = b->sift[i];
	//return distppf(a->sift, b->sift, SIFT_LENGTH);
	if (DIST_DESCS_EMV() > 0.5)
		return emvdistppf(a->sift, b->sift, SIFT_LENGTH);
	else if (DIST_DESCS_LP() != 2)
		return distlpf(a->sift, b->sift, SIFT_LENGTH, DIST_DESCS_LP());
	else
		return distppf(a->sift, b->sift, SIFT_LENGTH);
}

static double dist_descst(struct sift_keypoint *a, struct sift_keypoint *b,
		float t)
{
	if (DIST_DESCS_EMV() > 0.5) {
		float r = emvdistppf(a->sift, b->sift, SIFT_LENGTH);
		return r < t ? r : t;
	} else if (DIST_DESCS_LP() != 2) {
		float r=distlpf(a->sift, b->sift, SIFT_LENGTH,DIST_DESCS_LP());
		return r < t ? r : t;
	} else
		return distppft(a->sift, b->sift, SIFT_LENGTH, t);
}

static float euclidean_distance_topped(float *a, float *b, int n, float top)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		if (r > top)
			return top + 1;
		else
			r = hypot(r, b[i] - a[i]);
	return r;
}

static double sift_distance_topped(
		struct sift_keypoint *a,
		struct sift_keypoint *b,
		float top)
{
	return euclidean_distance_topped(a->sift, b->sift, SIFT_LENGTH, top);
}

static int nearestone(struct sift_keypoint *q, struct sift_keypoint *t, int nt, double *od)
{
	int bi = -1;
	double bd = INFINITY;
	FORI(nt) {
		double nb = dist_descs(q, t+i);
		if (nb < bd) {
			bd = nb;
			bi = i;
		}
	}
	assert(bi >= 0);
	*od = bd;
	return bi;
}

// returns i such that t[i] is as close as possible to q
static int fancynearest(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		double *od, double *opd)
{
	int besti = -1;
	double bestd = INFINITY, bestprev = bestd;
	FORI(nt) {
		double nb = dist_descs(q, t+i);
		if (nb < bestd) {
			bestprev = bestd;
			bestd = nb;
			besti = i;
		}
	}
	assert(besti >= 0);
	*od = bestd;
	*opd = bestprev;
	return besti;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int fancynearestt(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax)
{
	int besti = -1;
	float bestd = INFINITY;
	FORI(nt) {
		// TODO: be bold and change "dmax" to "bestd" in the call below
		float nb = dist_descst(q, t+i, dmax);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int fancynearestt_rad(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax, float dx, float dy)
{
	int besti = -1;
	float bestd = INFINITY;
	FORI(nt) {
		if (fabs(q->pos[0] - t[i].pos[0]) > dx) continue;
		if (fabs(q->pos[1] - t[i].pos[1]) > dy) continue;
		float nb = dist_descst(q, t+i, dmax);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int find_closest_keypoint(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax, float dx, float dy)
{
	int besti = -1;
	float bestd = dmax;
	FORI(nt) {
		if (fabs(q->pos[0] - t[i].pos[0]) > dx) continue;
		if (fabs(q->pos[1] - t[i].pos[1]) > dy) continue;
		float nb = sift_distance_topped(q, t+i, bestd);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

// returns i such that t[idxt[i]] is as close as possible to q
// if the distance is largest than dmax, return -1
static int find_closest_keypoint_idxs(struct sift_keypoint *q,
		struct sift_keypoint *t, int *idxt, int nit,
		float *od, float dmax, float dx, float dy)
{
	int besti = -1;
	float bestd = dmax;
	for (int i = 0; i < nit; i++)
	{
		int ii = idxt[i];
		if (fabs(q->pos[0] - t[ii].pos[0]) > dx) continue;
		if (fabs(q->pos[1] - t[ii].pos[1]) > dy) continue;
		float nb = sift_distance_topped(q, t+ii, bestd);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

// returns i such that t[idxt[i]] is as close as possible to q
static int fancynearest_idx(struct sift_keypoint *q,
		struct sift_keypoint *t, int *idxt, int nit,
		float *od, float *opd, float dmax)
{
	int besti = -1;
	float bestd = INFINITY, bestprev = bestd;
	for (int i = 0; i < nit; i++)
	{
		int ii = idxt[i];
		float nb = dist_descs(q, t+ii);
		if (nb < bestd) {
			bestprev = bestd;
			bestd = nb;
			besti = i;
		}
	}
	//assert(besti >= 0);
	*od = bestd;
	*opd = bestprev;
	return besti;
}

static
int (*siftlike_getpairs(
		struct sift_keypoint *ka, int na, 
		struct sift_keypoint *kb, int nb,
		int *np
		))[3]
{
	int (*p)[3] = xmalloc(na * sizeof * p);
	FORI(na) {
		double d;
		p[i][0] = i;
		p[i][1] = nearestone(ka+i, kb, nb, &d);
		p[i][2] = 10*log(d);
	}
	*np = na;
	return p;
}

static int compare_annpairs(const void *a, const void *b)
{
	float x = ((const struct ann_pair *)a)->v[0];
	float y = ((const struct ann_pair *)b)->v[0];
	int r = (x > y) - (x < y);
	//fprintf(stderr, "CMP %g %g = %d\n", x, y, r);
	return r;
}
static void sort_annpairs(struct ann_pair *t, int n)
{
	qsort(t, n, sizeof*t, compare_annpairs);
}

// get two lists of points, and produce a list of pairs
// (nearest match from a to b)
static
struct ann_pair *siftlike_get_annpairs(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *np
		)
{
	struct ann_pair *p = xmalloc(na * sizeof * p);
	FORI(na) {
		double d, dp;
		p[i].from = i;
		p[i].to = fancynearest(ka+i, kb, nb, &d, &dp);
		p[i].v[0] = d;
		p[i].v[1] = dp;
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, na);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*np = na;
	return p;
}

// get two lists of points, and produce a list of pairs
// (nearest match from a to b)
static
struct ann_pair *siftlike_get_annpairs_lowe(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *np,
		float loweratio
		)
{
	assert(loweratio < 1);
	assert(loweratio > 0);
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int cx = 0;
	FORI(na) {
		double d, dp;
		fancynearest(ka+i, kb, nb, &d, &dp);
		assert(dp >= d);
		if (d / dp < loweratio) {
			p[cx].from = i;
			p[cx].to = fancynearest(ka+i, kb, nb, &d, &dp);
			p[cx].v[0] = d;
			p[cx].v[1] = dp;
			cx += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, cx);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*np = cx;
	return p;
}

// get two lists of points, and produce a list of pairs
// (nearest match from a to b)
static
struct ann_pair *siftlike_get_annpairs_lowe2(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *np,
		float loweratio, float tdown, float tup
		)
{
	assert(loweratio < 1);
	assert(loweratio > 0);
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int cx = 0;
	int count_tdown = 0;
	int count_ratio = 0;
	int count_tup = 0;
	int count_ratiotup = 0;
	FORI(na) {
		double d, dp;
		fancynearest(ka+i, kb, nb, &d, &dp);
		assert(dp >= d);
		if ((d / dp < loweratio && d < tup) || d < tdown ) {
			p[cx].from = i;
			p[cx].to = fancynearest(ka+i, kb, nb, &d, &dp);
			p[cx].v[0] = d;
			p[cx].v[1] = dp;
			cx += 1;
		}
		if (d < tdown) count_tdown += 1;
		if (d/dp < loweratio) count_ratio += 1;
		if (d < tup) count_tup += 1;
		if (d/dp < loweratio && d < tup) count_ratiotup += 1;
	}
	fprintf(stderr, "count_tdown = %d\n", count_tdown);
	fprintf(stderr, "count_ratio = %d\n", count_ratio);
	fprintf(stderr, "count_tup = %d\n", count_tup);
	fprintf(stderr, "count_ratiotup = %d\n", count_ratiotup);
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, cx);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*np = cx;
	return p;
}

// get two lists of points, and produce a list of pairs
// (first nearest matches)
static
struct ann_pair *siftlike_get_accpairs(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *onp,
		float t
		)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float d;
		int j = fancynearestt(ka + i, kb, nb, &d, t);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

// get two lists of points, and produce a list of pairs
// (first nearest matches)
static
struct ann_pair *siftlike_get_accpairsrad(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *onp,
		float t, float dx, float dy
		)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float d;
		int j = fancynearestt_rad(ka + i, kb, nb, &d, t, dx, dy);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

//#define VERBOSE true
#include "ok_list.c"
#include "grid.c"

struct ok_grid {
	struct ok_list l[1];
	struct grid g[1];
	int *buf;
};

static void ok_grid_init(struct ok_grid *o, int np,
		float x0[2], float dx[2], int n[2])
{
	int nr = n[0] * n[1];
	ok_init(o->l, nr, np);
	grid_init(o->g, 2, x0, dx, n);
	o->buf = xmalloc(np*sizeof*o->buf);
}

static void ok_grid_free(struct ok_grid *o)
{
	ok_free(o->l);
	free(o->buf);
}

static int ok_grid_add_point(struct ok_grid *o, int p, float x[2])
{
	int r = grid_locate(o->g, x);
	ok_add_point(o->l, r, p);
	return r;
}

static int ok_neighboring_points(struct ok_grid *o, float x[2])
{
	int r[4], nr = grid_locate_overlapping(r, o->g, x);
	assert(nr <= 4);
	//for (int i=0;i<nr;i++) fprintf(stderr, "\trglos{%g %g} [%d:%d] = %d\n", x[0], x[1], i, nr, r[i]);
	int cx = 0;
	for (int i = 0; i < nr; i++)
	{
		int nri = ok_which_points(o->l, r[i]);
		for (int j = 0; j < nri; j++)
			o->buf[cx++] = o->l->buf[j];
	}
	return cx;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int find_closest_in_grid(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax, struct ok_grid *g)
{
	// build the list of neighbors to traverse
	int nn = ok_neighboring_points(g, q->pos);
	int *nbuf = g->buf;
	//fprintf(stderr, "site ( %g , %g ) has %d neighbors\n", q->pos[0], q->pos[1], nn);
	//for(int i = 0; i < nn; i++) fprintf(stderr, "\t%d\n", nbuf[i]);
	
	// compute the closest point on this list
	int besti = -1;
	float bestd = dmax;
	FORI(nn) {
		int idx = nbuf[i];
		if (fabs(q->pos[0] - t[idx].pos[0]) > g->g->dx[0]) continue;
		if (fabs(q->pos[1] - t[idx].pos[1]) > g->g->dx[1]) continue;
		float nb = sift_distance_topped(q, t+idx, bestd);
		if (nb < bestd) {
			bestd = nb;
			besti = idx;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

static
struct ann_pair *compute_sift_matches_locally(int *onp,
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		float t, float dx, float dy, int w, int h)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }

	// build a grid structure for each list of keypoints
	// NOTE: only the grid of "ka" is actually used
	float x0[2] = {0, 0};
	float dxy[2] = {dx, dy};
	int n[2] = {1+(w-1)/dx, 1+(h-1)/dy};
	//int nr = n[0] * n[1];
	//struct ok_grid ga[1]; ok_grid_init(ga, na, x0, dxy, n);
	struct ok_grid gb[1]; ok_grid_init(gb, nb, x0, dxy, n);
	//for (int i = 0; i < na; i++) ok_grid_add_point(ga, i, ka[i].pos);
	for (int i = 0; i < nb; i++) ok_grid_add_point(gb, i, kb[i].pos);
	//ok_display_tables(gb->l);
	//ok_assert_consistency(gb->l);

	// compute the pairs
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	for (int i = 0; i < na; i++) {
		float d;
		int j = find_closest_in_grid(ka + i, kb, nb, &d, t, gb);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	sort_annpairs(p, np);

	// free the grid structures
	//ok_grid_free(ga);
	ok_grid_free(gb);

	// return values
	*onp = np;
	return p;
}

static
struct ann_pair *compute_sift_matches(int *onp,
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		float t, float dx, float dy)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	for (int i = 0; i < na; i++) {
		float d;
		int j = find_closest_keypoint(ka + i, kb, nb, &d, t, dx, dy);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	sort_annpairs(p, np);
	*onp = np;
	return p;
}

static void get_bbx(double out_min[2], double out_max[2],
		struct sift_keypoint *p, int np)
{
	out_min[0] = out_min[1] = INFINITY;
	out_max[0] = out_max[1] = -INFINITY;
	for (int i = 0; i < np; i++)
	for (int l = 0; l < 2; l++)
	{
		out_min[l] = fmin(out_min[l], p[i].pos[l]);
		out_max[l] = fmax(out_max[l], p[i].pos[l]);
	}
}


static void vector_times_matrix(double xA[3], double x[3], double A[9])
{
	xA[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
	xA[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
	xA[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}

// compute the scalar product of two vectors
static double scalar_product(double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// cut a line with a segment (returns true if they cut)
static bool cut_line_with_segment(double out[2], double line[3],
		double p[2], double q[2])
{
	// points in "oriented" projective coordinates
	double pp[3] = {p[0], p[1], 1};
	double qq[3] = {q[0], q[1], 1};

	// sign of each point (says on which side of the line each point is)
	double sp = scalar_product(pp, line);
	double sq = scalar_product(qq, line);

	// if signs are different, the line crosses the segment
	if (sp * sq < 0) {
		// line trough points p and q
		double pq[3]; vector_product(pq, pp, qq);

		// intersection of "line" and "pq"
		double ii[3]; vector_product(ii, pq, line);

		// recover affine coordinates
		out[0] = ii[0] / ii[2];
		out[1] = ii[1] / ii[2];
		return true;
	}
	return false;
}

// cut a line with a rectangle (returns true if they cut)
static bool cut_line_with_rectangle(double out_a[2], double out_b[2],
		double line[3], double rec_from[2], double rec_to[2])
{
	//double nnn = hypot(line[0], line[1]);
	//fprintf(stderr, "clwr (%g %g %g) (%g %g)-(%g %g)...",
	//		line[0]/nnn, line[1]/nnn, line[2]/nnn,
	//		rec_from[0], rec_from[1], rec_to[0], rec_to[1]);
	// four vertices of the rectangle
	double v[4][2] = {
		{ rec_from[0], rec_from[1] },
		{ rec_to[0]  , rec_from[1] },
		{ rec_to[0]  , rec_to[1]   },
		{ rec_from[0], rec_to[1]   }
	};

	// intersections with each of the edges
	bool xP[4]; // whether it intersects
	double x[4][2]; // where it intersects
	for (int i = 0; i < 4; i++)
		xP[i] = cut_line_with_segment(x[i], line, v[i], v[ (i+1)%4 ] );

	// write output
	int n_intersections = xP[0] + xP[1] + xP[2] + xP[3];
	if (n_intersections == 2) { // generic case: 2 intersections
		int cx = 0;
		for (int i = 0; i < 4; i++)
			if (xP[i])
			{
				double *out = cx ? out_b : out_a;
				out[0] = x[i][0];
				out[1] = x[i][1];
				cx += 1;
			}
		//fprintf(stderr, "cx=%d out=(%g %g)-(%g %g)\n",cx, out_a[0], out_a[1], out_b[0], out_b[1]);
		return true;
	}
	//fprintf(stderr, "nothing\n");
	return false;
}

// draw a segment between two points
static void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px+j*slope), j+py, e);
		}
	}
}

static
void traverse_segment_thick_precise(
		double px, double py, double qx, double qy,
		void (*f)(int,int,void*), void *e)
{
	if (qx + qy < px + py) // bad quadrants
		traverse_segment_thick_precise(qx, qy, px, py, f, e);
	else {
		if (fabs(qx - px) > qy - py) { // horitzontal
			double slope = (qy - py); slope /= (qx - px);
			assert(px <= qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i <= qx-px; i++) {
				double exact = py + i*slope;
				int whole = lrint(exact);
				//double part = fabs(whole - exact);
				//int owhole = (whole<exact)?whole+1:whole-1;
				//assert(part <= 0.5);
				f(i+px, whole, e);
				f(i+px, whole+1, e);
				f(i+px, whole-1, e);
				f(i+px, whole+2, e);
				f(i+px, whole-2, e);
				//f(i+px, whole, 1-part, e);
				//f(i+px, owhole, part, e);
			}
		} else { // vertical
			double slope = (qx - px); slope /= (qy - py);
			assert(fabs(qy - py) >= fabs(qx - px));
			assert(py <= qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++) {
				double exact = px + j*slope;
				int whole = lrint(exact);
				//double part = fabs(whole - exact);
				//int owhole = (whole<exact)?whole+1:whole-1;
				//assert(part <= 0.5);
				f(whole, j+py, e);
				f(whole+1, j+py, e);
				f(whole-1, j+py, e);
				f(whole+2, j+py, e);
				f(whole-2, j+py, e);
				//f(whole, j+py, 1-part, e);
				//f(owhole, j+py, part, e);
			}
		}
	}
}

// draw a segment between two points (somewhat anti-aliased)
static
void traverse_segment_thick(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, qx, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment_thick(qx, qy, px, py, f, e);
	else {
		//fprintf(stderr, "tsthick (%d %d)-(%d %d)\n", px, py, qx, qy);
		if (abs(qx - px) > qy - py) { // horitzontal
			float slope = (qy - py); slope /= (qx - px);
			assert(px <= qx);
			assert(fabs(slope) <= 1);
			for (int i = 0; i <= qx-px; i++) {
				float exact = py + i*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				//int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(i+px, whole, e);
				f(i+px, whole+1, e);
				f(i+px, whole-1, e);
				//f(i+px, whole, 1-part, e);
				//f(i+px, owhole, part, e);
			}
		} else { // vertical
			float slope = (qx - px); slope /= (qy - py);
			assert(abs(qy - py) >= abs(qx - px));
			assert(py <= qy);
			assert(fabs(slope) <= 1);
			for (int j = 0; j <= qy-py; j++) {
				float exact = px + j*slope;
				int whole = lrint(exact);
				float part = fabs(whole - exact);
				//int owhole = (whole<exact)?whole+1:whole-1;
				assert(part <= 0.5);
				f(whole, j+py, e);
				f(whole+1, j+py, e);
				f(whole-1, j+py, e);
				//f(whole, j+py, 1-part, e);
				//f(owhole, j+py, part, e);
			}
		}
	}
}

struct grille_traversal_state {
	int *buf, nbuf, maxbuf;
	int grille_width, grille_height;
	struct ok_grid *o;
	float eps;
	float img_line[3];
	struct sift_keypoint *kb;
};

static
float signed_distance_point_to_line(float l[3], float x[2])
{
	return (l[0] * x[0] + l[1] * x[1] + l[2]) / hypot(l[0], l[1]);
}

// add to the buffer the indexes of keypoints found in grille position (i,j)
static
void grille_traversal_function(int i, int j, void *ee)
{
	struct grille_traversal_state *e = ee;
	//fprintf(stderr, "GTF CALL\n");fflush(stderr);
	//fprintf(stderr, "grille pos %d %d (%g)\n", i, j, a);

	//if (i < 0) i = 0;
	//if (j < 0) j = 0;
	//if (i >= e->grille_width ) i = e->grille_width  - 1;
	//if (j >= e->grille_height) j = e->grille_height - 1;
	if (i < 0 || j < 0) return;
	if (i >= e->grille_width ) return;
	if (j >= e->grille_height) return;
	int ridx = j * e->grille_width + i;
	int np = ok_which_points(e->o->l, ridx);
	//fprintf(stderr, "\t%d points here\n", np);
	if (e->nbuf + np >= e->maxbuf)
		fail("grille buffer overflow!");
	for (int k = 0; k < np; k++)
	{
		int idx = e->o->l->buf[k];
		float d = signed_distance_point_to_line(
				e->img_line, e->kb[idx].pos);
		if (fabs(d) < e->eps)
			e->buf[e->nbuf++] = idx;
		//e->buf[e->nbuf++] = e->o->l->buf[k];
	}
}

#include "iio.h"

static int insideP(int w, int h, int i, int j)
{
	return i>=0 && j>= 0 && i<w && j<h;
}

struct plot_state {
	uint8_t *x;
	int w, h, r, g, b;
};

static void pixel_plotter(int i, int j, void *ee)
{
	struct plot_state *e = ee;
	if (insideP(e->w, e->h, i, j))
	{
		e->x[3*(e->w*j+i)+0] = e->r;
		e->x[3*(e->w*j+i)+1] = e->g;
		e->x[3*(e->w*j+i)+2] = e->b;
	}
}

static
void overlay_red_thick_line(uint8_t *x, int w, int h, double a[2], double b[2])
{
	struct plot_state e[1];
	e->x = x;
	e->w = w;
	e->h = h;
	e->r = 255;
	e->g = 0;
	e->b = 0;
	traverse_segment_thick(
			round(a[0]), round(a[1]),
			round(b[0]), round(b[1]),
			pixel_plotter, e);
}

static
struct ann_pair *sift_fm_pairs(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		float tau, double fm[9], float eps,
		int *out_np)
{
	if (na == 0 || nb == 0) { *out_np=0; return NULL; }
	struct ann_pair *p = xmalloc((na*nb) * sizeof * p);
	int num_pairs = 0;

	// compute bounding box of points on image B
	double minxy[2], maxxy[2];
	get_bbx(minxy, maxxy, kb, nb);
	//fprintf(stderr, "B's bbx = (%g %g)-(%g %g)\n", minxy[0], minxy[1], maxxy[0], maxxy[1]);

	// prepare occupancy list of image B
	float dxy[2] = {eps, eps};
	int nxy[2] = {1+(maxxy[0]-minxy[0])/eps, 1+(maxxy[1]-minxy[1])/eps};
	struct ok_grid gb[1];
	ok_grid_init(gb, nb,(float[]){minxy[0],minxy[1]}, dxy, nxy);
	for (int i = 0; i < nb; i++)
		ok_grid_add_point(gb, i, kb[i].pos);

	// transformation from "quadrille" to "pixel" coordinates
	int qw = nxy[0];
	int qh = nxy[1];
	int q[qh][qw];
	for (int j = 0; j < qh; j++)
	for (int i = 0; i < qw; i++)
		q[j][i] = j * qw + i;
	//double aff_a[2] = { 1/eps, 1/eps};
	//double aff_b[2] = { -minxy[0]/eps, -minxy[1]/eps};

	// buffer for point indexes
	int buf[nb];
	struct grille_traversal_state e[1];
	e->buf = buf;
	e->nbuf = 0;
	e->maxbuf = nb;
	e->o = gb;
	e->grille_width = qw;
	e->grille_height = qh;
	e->eps = eps/2;
	e->kb = kb;

	fprintf(stderr, "quadrille (%d %d) with eps=%g\n", qw, qh, eps);

	// for each A point, etc
	//assert(aff_a[0] == aff_a[1]);
	for (int i = 0; i < na; i++)
	{
		e->nbuf = 0;
		double x[3] = { ka[i].pos[0], ka[i].pos[1], 1};
		//fprintf(stderr, "treating %dth A point (%g %g)\n",i,x[0],x[1]);
		double epix[3]; // epipolar line in right IMAGE coordinates
		vector_times_matrix(epix, x, fm);
		double factor = hypot(epix[0], epix[1]);
		for (int l = 0; l < 3; l++) epix[l] /= factor;
		for (int l = 0; l < 3; l++) e->img_line[l] = epix[l];
		double epiX[3]; // epipolar line in right GRILLE coordinates
		epiX[0] = epix[0];
		epiX[1] = epix[1];
		epiX[2] = (epix[2] + epix[0]*minxy[0] + epix[1]*minxy[1]) / eps;
		//epiX[2] = (epix[2] + aff_b[0]*epix[0] + aff_b[1]*epix[1]) * aff_a[0];
		//epiX[2] = (epix[2] + aff_b[0]*epix[0] + aff_b[1]*epix[1]) * aff_a[0];
		//fprintf(stderr, "epix = %g %g %g\n", epix[0], epix[1], epix[2]);
		//fprintf(stderr, "epiX = %g %g %g\n", epiX[0], epiX[1], epiX[2]);
		double gfrom[2] = {0, 0};
		double gto[2] = {qw, qh};
		double pa[2], pb[2];
		cut_line_with_rectangle(pa, pb, epix, minxy, maxxy);
		if (cut_line_with_rectangle(pa, pb, epiX, gfrom, gto))
		{
			//int ipa[2] = {round(pa[0]), round(pa[1])};
			//int ipb[2] = {round(pb[0]), round(pb[1])};
			traverse_segment_thick_precise(
					pa[0], pa[1], pb[0], pb[1],
					grille_traversal_function, e);
		}
		//fprintf(stderr, "enbuf = %d\n", e->nbuf);
		// now "e->buf" contains the list
		// of the "e->nbuf" candidate keypoints (from image B)
		// these keypoints must be compared to the keypoint "ka[i]"
		// and the best one, if any, added to the list "p"
		//
		if (isnan(tau) && i == 0) { // debug mode
			for (int k = 0; k < e->nbuf; k++)
			{
				int j = e->buf[k];
				assert(j >= 0);
				assert(j < nb);
				p[num_pairs].from = i;
				p[num_pairs].to = j;
				p[num_pairs].v[0] = 0;
				p[num_pairs].v[1] = NAN;
				//fprintf(stderr, "added %dth pair {%d} [%d %d] (d=%g)\n", num_pairs, cidx, i, j, dist);
				//if (dist_point_line(epix, kb[j].pos[0], kb[j].pos[1]) < eps)
					num_pairs += 1;
			}
			int dbg_w = maxxy[0];
			int dbg_h = maxxy[1];
			uint8_t *dbg = xmalloc(3 * dbg_w * dbg_h);
			for (int k = 0; k < dbg_w * dbg_h * 3; k++)
				dbg[k] = 0;
			cut_line_with_rectangle(pa, pb, epix, minxy, maxxy);
			fprintf(stderr, "pa = %g %g\n", pa[0], pa[1]);
			fprintf(stderr, "pb = %g %g\n", pb[0], pb[1]);
			overlay_red_thick_line(dbg, dbg_w, dbg_h, pa, pb);
			for (int k = 0; k < num_pairs; k++)
			{
				int idx = p[k].to;
				int xp = round(kb[idx].pos[0]);
				int yp = round(kb[idx].pos[1]);
				if (insideP(dbg_w, dbg_h, xp, yp))
					dbg[3*(dbg_w*yp+xp)+1] = 255;
			}
			iio_write_image_uint8_vec("/tmp/dbg_sfm.png", dbg, dbg_w, dbg_h, 3);
			free(dbg);
			goto acabemaqui;
		}
		float dist, dsecond;
		int cidx;//
		if (tau < 0 && tau > -1) {
			cidx = fancynearest_idx(ka + i,
					kb, e->buf, e->nbuf,
					&dist, &dsecond, 500);
			if (dist > 500 || (dist / dsecond > -tau))
				cidx = -1;
		}
		else
			cidx = find_closest_keypoint_idxs(ka + i,
				kb, e->buf, e->nbuf,
				&dist, tau, INFINITY, INFINITY);
		if (cidx >= 0)
		{
			int j = e->buf[cidx];
			assert(j >= 0);
			assert(j < nb);
			p[num_pairs].from = i;
			p[num_pairs].to = j;
			p[num_pairs].v[0] = dist;
			p[num_pairs].v[1] = NAN;
			//fprintf(stderr, "added %dth pair {%d} [%d %d] (d=%g)\n", num_pairs, cidx, i, j, dist);
			num_pairs += 1;
		}
	}

acabemaqui:
	sort_annpairs(p, num_pairs);
	*out_np = num_pairs;
	return p;
}


// get three lists of points, and produce a list of matching triplets
// (first nearest matches)
static
struct ann_trip *siftlike_get_triplets(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		struct sift_keypoint *kc, int nc,
		int *onp,
		float t
		)
{
	if (na == 0 || nb == 0 || nc == 0) { *onp=0; return NULL; }
	struct ann_trip *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float db, dc;
		int jb = fancynearestt(ka + i, kb, nb, &db, t);
		int jc = fancynearestt(ka + i, kc, nc, &dc, t);
		if (jb >= 0 && jc >= 0) {
			p[np].froma = i;
			p[np].tob = jb;
			p[np].toc = jc;
			p[np].v[0] = hypot(db,dc);
			p[np].v[1] = db;
			p[np].v[2] = dc;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	//sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

// get three lists of points, and produce a list of matching triplets
// (first nearest matches)
static
struct ann_trip *siftlike_get_tripletsrad(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		struct sift_keypoint *kc, int nc,
		int *onp,
		float t, float rx, float ry
		)
{
	if (na == 0 || nb == 0 || nc == 0) { *onp=0; return NULL; }
	struct ann_trip *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float db, dc;
		int jb = fancynearestt_rad(ka + i, kb, nb, &db, t, rx, ry);
		int jc = fancynearestt_rad(ka + i, kc, nc, &dc, t, rx, ry);
		if (jb >= 0 && jc >= 0) {
			p[np].froma = i;
			p[np].tob = jb;
			p[np].toc = jc;
			p[np].v[0] = hypot(db,dc);
			p[np].v[1] = db;
			p[np].v[2] = dc;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	//sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

// get two lists of points, and produce a list of pairs
// (all nearest matches)
static
struct ann_pair *siftlike_get_allpairs(
		struct sift_keypoint *ka, int na, 
		struct sift_keypoint *kb, int nb,
		int *onp,
		float t
		)
{
	struct ann_pair *p = NULL;//xmalloc(na * sizeof * p);
	int np = 0, np_top = 0;
	FORI(na) FORJ(nb) {
		//double d = dist_descs(ka+i, kb+j);
		double d = dist_descst(ka+i, kb+j, t);
		if (d < t) {
			if (np >= np_top) {
				np_top = 0x100 + 2 * (np_top + 1);
				p = xrealloc(p, np_top * sizeof * p);
			}
			struct ann_pair *pi = p + np;
			pi->from = i;
			pi->to = j;
			pi->v[0] = d;
			pi->v[1] = NAN;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	fprintf(stderr, "NOW, to sort the %d pairs\n", np);
	sort_annpairs(p, np);
	FILE *g = xfopen("/tmp/siftpairsd.txt", "w");
	FORI(np) fprintf(g, "%g\n", p[i].v[0]);
	xfclose(g);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

#if USE_KDTREE

#define POINTY_HACK 4387
static void *int_to_pointer(int x)
{
	void *p = (void *)(x + POINTY_HACK);
	return p;
}
static int pointer_to_int(void *p)
{
	int x = -POINTY_HACK + (int)p;
	return x;
}

// using kdtrees,
// get two lists of points, and produce a list of pairs
// (nearest match from a to b)
#endif
static
struct ann_pair *siftlike_get_annpairs_kd(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *np
		)
{
#if USE_KDTREE
	struct ann_pair *p = xmalloc(na * sizeof * p);
	struct kdtree *kd = kd_create(SIFT_LENGTH);
	fprintf(stderr, "inserting points on kdtree\n");
	FORI(nb)
		kd_insertf(kd, kb[i].sift, int_to_pointer(i));
	fprintf(stderr, "\t...done\n");
	fprintf(stderr, "finding nearest neighbors\n");
	FORI(na) {
		float *f = ka[i].sift;
		struct kdres *res = kd_nearestf(kd, f);
		if (1 == kd_res_size(res)) {
			float fn[SIFT_LENGTH];
			void *pi = kd_res_itemf(res, fn);
			if (!pi) fail("epa aquí!");
			p[i].from = i;
			p[i].to = pointer_to_int(pi);
			p[i].v[0] = distppf(fn, f, SIFT_LENGTH);
			p[i].v[1] = INFINITY;
		} else
			fail("no en tenim cap");
	}
	fprintf(stderr, "\t...done\n");
	kd_free(kd);

	sort_annpairs(p, na);
	*np = na;
	return p;
#else
	return siftlike_get_annpairs(ka, na, kb, nb, np);
#endif
}

//int compare_annpairs(const struct ann_pair *a, const struct ann_pair *b)

SMART_PARAMETER(PSIFT_DISTTHRE,-1)
SMART_PARAMETER(PSIFT_LOWERATIO,-1)

// mask pairs by using some thresholds
static
void siftlike_maskpairs(bool *mask, struct ann_pair *t, int n)
{
	FORI(n) mask[i] = true;
	float dthre = PSIFT_DISTTHRE();
	float lrat = PSIFT_LOWERATIO();
	if (dthre > 0)
		FORI(n)
			if (t[i].v[0] > dthre)
				mask[i] = false;
	if (lrat > 0)
		FORI(n)
			if (t[i].v[0] > lrat * t[i].v[1])
				mask[i] = false;
}


static float sift_d2(struct sift_keypoint *a, struct sift_keypoint *b)
{
	float dx = b->pos[0] - a->pos[0];
	float dy = b->pos[1] - a->pos[1];
	return hypot(dx,dy);
}

static float sift_d128(struct sift_keypoint *a, struct sift_keypoint *b)
{
	return dist_descs(a, b);
	//return distppf(a->sift, b->sift, SIFT_LENGTH);
}

// TODO: make this function run in linear or near-linear time
// (by assuming that the points are uniformly distributed on the plane, and
// using an efficient data structure for localizing them)
// Right now, it is prohibitively slow for the common case of more than 10.000
// points
static
void sift_remove_redundancy(bool *mask,
		struct sift_keypoint *t, int n,
		float dist_plane, float dist_sift)
{
	FORI(n) mask[i] = true;

	FORI(n) FORJ(i)
		if (mask[i])
			if (sift_d2(t+i, t+j) < dist_plane)
				if (sift_d128(t+i, t+j) < dist_sift)
					mask[j] = false;
}

static
int splitloc(int *t, int w, float rx, float ry, float ox, float oy, float *x)
{
	int r = 0;
	int ix = (int)(x[0]/(rx-ox));
	int iy = (int)(x[1]/(ry-oy));
	fprintf(stderr, "\tix=%d, rx=%g, ox=%g, x[0]=%g\n", ix,rx,ox,x[0]);
	fprintf(stderr, "\tiy=%d, ry=%g, oy=%g, x[1]=%g\n", iy,ry,oy,x[1]);
	assert(ix >= 0);
	assert(ix < w);
	assert(iy >= 0);
	assert(ix * (rx-ox) <= x[0]);
	assert(iy * (ry-oy) <= x[1]);
	assert((ix+1) * (rx-ox) >= x[0]);
	assert((iy+1) * (ry-oy) >= x[1]);
	t[r++] = iy*w + ix;
	bool ovx = ix > 0 && ix * (rx - ox) + ox >= x[0];
	bool ovy = iy > 0 && iy * (ry - oy) + oy >= x[1];
	if (ovx) t[r++] = iy*w + ix - 1;
	if (ovy) t[r++] = (iy-1)*w + ix;
	if (ovx && ovy) t[r++] = (iy-1)*w + ix - 1;
	fprintf(stderr, "\t\tt={%d %d %d %d}\n", t[0], t[1], t[2], t[3]);
	return r;
}

static
int siftsplit(struct sift_keypoint *p, int n,
		float rx, float ry, float ox, float oy,
		int (*mask)[5])
{
	FORI(n) mask[i][0] = 0;
	FORI(n) FORJ(4) mask[i][1+j] = -1;
	float maxx = -INFINITY; FORI(n) maxx = fmax(maxx,p[i].pos[0]);
	float maxy = -INFINITY; FORI(n) maxy = fmax(maxy,p[i].pos[1]);
	//int w = calcula_amplada(...);
	//int r = calcula_nombre_de_rectangles(rx,ry,ox,oy,maxx,maxy);
	int w = 1+(int)(maxx/(rx-ox));
	int r = (1+w) * (1+(int)(maxy/(ry-oy)));
	fprintf(stderr, "w = %d, r = %d\n", w, r);
	FORI(n) {
		fprintf(stderr, "p[%d] = %g , %g\n",i,p[i].pos[0],p[i].pos[1]);
		mask[i][0] = splitloc(mask[i]+1, w, rx, ry, ox, oy, p[i].pos);
		fprintf(stderr, "belongs to %d region:\n", mask[i][0]);
		FORJ(mask[i][0])
			fprintf(stderr, "\t%d\n", mask[i][1+j]);
		fprintf(stderr, "\n");
	}
	return r;
}

static void affine_mapf(float y[2], float A[6], float x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

static void homographic_mapf(float y[2], float H[9], float x[2])
{
	float z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	y[1] = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
	//y[0] = x[0];
	//y[1] = x[1];
}

static
void siftaff(struct sift_keypoint *t, int n, float A[9])
{
	float det = A[0]*A[4] - A[1]*A[3];
	//fprintf(stderr, "det = %g\n", det);
	FORI(n) {
		struct sift_keypoint *k = t+i;
		float vec[2] = {cos(k->orientation), sin(k->orientation)};
		float x[2], rvec[2];
		affine_mapf(x, A, k->pos);
		rvec[0] = A[0]*vec[0] + A[1]*vec[1];
		rvec[1] = A[3]*vec[0] + A[4]*vec[1];
		FORJ(2) k->pos[j] = x[j];
		k->scale *= det;
		k->orientation = atan2(rvec[1], rvec[0]);
	}
}

static
void sifthom(struct sift_keypoint *t, int n, float H[9])
{
	// TODO XXX ERROR FIXME : update the scale and orientation accordingly!
	FORI(n) {
		struct sift_keypoint *k = t+i;
		float x[2];
		homographic_mapf(x, H, k->pos);
		FORJ(2) k->pos[j] = x[j];
	}
}
