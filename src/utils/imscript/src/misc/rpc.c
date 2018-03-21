// rational polynomial coefficient stuff

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "xfopen.c"



// rational polynomial coefficients (and those of the inverse model)
struct rpc {
	double numx[20];
	double denx[20];
	double numy[20];
	double deny[20];
	double scale[3], offset[3];

	double inumx[20];
	double idenx[20];
	double inumy[20];
	double ideny[20];
	double iscale[3], ioffset[3];

	double dmval[4];
	double imval[4];
};

// set all the values of an rpc model to NAN
static void nan_rpc(struct rpc *p)
{
	int nd = sizeof*p/sizeof(double);
	double *t = (double*)p;
	for (int i = 0; i < nd; i++)
		t[i] = NAN;
}

// like strcmp, but finds a needle
static int strhas(char *haystack, char *needle)
{
	char *r = strstr(haystack, needle);
	return r ? 0 : 1;
}

// returns the first character not in the prefix, or NULL if no prefix
static char *prefix(char *s, char *p)
{
	int i = 0;
	while (s[i] && p[i] && s[i] == p[i])
		i++;
	if (!s[i] || p[i])
		return NULL;
	else
		return s + i;
}

// if  s is the form "%s_%d"%(s,p), return %d bounded on [0,19]
static int pix(char *s, char *p)
{
	char *idx = prefix(s, p);
	if (!idx) return -1;
	int r = atoi(idx) - 1;
	if (r < 0) r = 0;
	if (r > 19) r = 19;
	return r;
}

// add a value to the rpc
static void add_tag_to_rpc(struct rpc *p, char *tag, double x)
{
	int t;
	if (false);
	else if (0 == strhas(tag, "SAMP_OFF"))          p->offset [0] = x;
	else if (0 == strhas(tag, "SAMP_SCALE"))        p->scale  [0] = x;
	else if (0 == strhas(tag, "LINE_OFF"))          p->offset [1] = x;
	else if (0 == strhas(tag, "LINE_SCALE"))        p->scale  [1] = x;
	else if (0 == strhas(tag, "HEIGHT_OFF"))        p->offset [2] = x;
	else if (0 == strhas(tag, "HEIGHT_SCALE"))      p->scale  [2] = x;
	else if (0 == strhas(tag, "LONG_OFF"))          p->ioffset[0] = x;
	else if (0 == strhas(tag, "LONG_SCALE"))        p->iscale [0] = x;
	else if (0 == strhas(tag, "LAT_OFF"))           p->ioffset[1] = x;
	else if (0 == strhas(tag, "LAT_SCALE"))         p->iscale [1] = x;
	else if (0 == strhas(tag, "HEIGHT_OFF"))        p->ioffset[2] = x;
	else if (0 == strhas(tag, "HEIGHT_SCALE"))      p->iscale [2] = x;
	else if (0 <= (t=pix(tag, "SAMP_NUM_COEFF_")))  p->numx   [t] = x;
	else if (0 <= (t=pix(tag, "SAMP_DEN_COEFF_")))  p->denx   [t] = x;
	else if (0 <= (t=pix(tag, "LINE_NUM_COEFF_")))  p->numy   [t] = x;
	else if (0 <= (t=pix(tag, "LINE_DEN_COEFF_")))  p->deny   [t] = x;
	else if (0 <= (t=pix(tag, "iSAMP_NUM_COEFF_"))) p->inumx  [t] = x;
	else if (0 <= (t=pix(tag, "iSAMP_DEN_COEFF_"))) p->idenx  [t] = x;
	else if (0 <= (t=pix(tag, "iLINE_NUM_COEFF_"))) p->inumy  [t] = x;
	else if (0 <= (t=pix(tag, "iLINE_DEN_COEFF_"))) p->ideny  [t] = x;
	else if (0 == strhas(tag, "FIRST_ROW"))         p->dmval  [1] = x;
	else if (0 == strhas(tag, "FIRST_COL"))         p->dmval  [0] = x;
	else if (0 == strhas(tag, "LAST_ROW"))          p->dmval  [3] = x;
	else if (0 == strhas(tag, "LAST_COL"))          p->dmval  [2] = x;
	else if (0 == strhas(tag, "FIRST_LON"))         p->imval  [0] = x;
	else if (0 == strhas(tag, "FIRST_LAT"))         p->imval  [1] = x;
	else if (0 == strhas(tag, "LAST_LON"))          p->imval  [2] = x;
	else if (0 == strhas(tag, "LAST_LAT"))          p->imval  [3] = x;
}

// process an input line
static double get_xml_tagged_number(char *tag, char *line)
{
	char buf[0x100], buf2[0x100];
	double x;
	int r = sscanf(line, " <%[^>]>%lf</%[^>]>", buf, &x, buf2);
	if (r == 3) {
		strcpy(tag, buf);
		return x;
	} else return NAN;
}

// read an XML file specifying an RPC model
void read_rpc_file_xml(struct rpc *p, char *filename)
{
	nan_rpc(p);
	FILE *f = xfopen(filename, "r");
	int n = 0x100, o = 1;
	while (1) {
		char line[n], tag[n], *sl = fgets(line, n, f);
		if (!sl) break;
		if (o && 0 == strhas(line, "<Inverse_Model>")) o = 0;
		tag[0] = 'i';
		double x = get_xml_tagged_number(tag+1, line);
		if (isfinite(x)) {
			//fprintf(stderr, "%s [%d]: %g\n", tag+o, o, x);
			add_tag_to_rpc(p, tag+o, x);
		}
	}
	xfclose(f);
	p->ioffset[2] = p->offset[2];
	p->iscale[2] = p->scale[2];
}

#define FORI(n) for (int i = 0; i < (n); i++)

void print_rpc(FILE *f, struct rpc *p, char *n)
{
	FORI(20) fprintf(f, "rpc(%s) numx[%d] = %.18lf\n",n,i,p->numx[i]);
	FORI(20) fprintf(f, "rpc(%s) denx[%d] = %.18lf\n",n,i,p->denx[i]);
	FORI(20) fprintf(f, "rpc(%s) numy[%d] = %.18lf\n",n,i,p->numy[i]);
	FORI(20) fprintf(f, "rpc(%s) deny[%d] = %.18lf\n",n,i,p->deny[i]);
	FORI(20) fprintf(f, "rpc(%s) inumx[%d] = %.18lf\n",n,i,p->inumx[i]);
	FORI(20) fprintf(f, "rpc(%s) idenx[%d] = %.18lf\n",n,i, p->idenx[i]);
	FORI(20) fprintf(f, "rpc(%s) inumy[%d] = %.18lf\n",n,i,p->inumy[i]);
	FORI(20) fprintf(f, "rpc(%s) ideny[%d] = %.18lf\n",n,i,p->ideny[i]);
	FORI(3)fprintf(stderr,"rpc(%s) scale[%d] = %.18lf\n",n,i,p->scale[i]);
	FORI(3)fprintf(stderr,"rpc(%s) offset[%d] = %.18lf\n",n,i,p->offset[i]);
	FORI(3)fprintf(stderr,"rpc(%s) iscale[%d] = %.18lf\n",n,i,p->iscale[i]);
	FORI(3)fprintf(stderr,"rpc(%s) ioffs[%d] = %.18lf\n",n,i,p->ioffset[i]);

}

// evaluate a polynomial of degree 3
double eval_pol20(double c[20], double x, double y, double z)
{
	// XXX WARNING: inversion here!
	//double col = y;
	//double lig = x;
	//double alt = z;
	//double m[20] = {1, lig, col, alt, lig*col,
	//	lig*alt, col*alt, lig*lig, col*col, alt*alt,
	//	col*lig*alt, lig*lig*lig, lig*col*col, lig*alt*alt, lig*lig*col,
	//	col*col*col, col*alt*alt, lig*lig*alt, col*col*alt, alt*alt*alt
	//};
	double m[20] = {1, x, y, z, x*y,
		x*z, y*z, x*x, y*y, z*z,
		x*y*z, x*x*x, x*y*y, x*z*z, x*x*y,
		y*y*y, y*z*z, x*x*z, y*y*z, z*z*z};
	double r = 0;
	for (int i = 0; i < 20; i++)
		r += c[i]*m[i];
	return r;
}

double eval_pol20_dx(double c[20], double x, double y, double z)
{
	double m[20] = {0, 1, 0, 0, y,
		z, 0, 2*x, 0, 0,
		y*z, 3*x*x, y*y, z*z, 2*x*y,
		0, 0, 2*x*z, 0, 0};
	double r = 0;
	for (int i = 0; i < 20; i++)
		r += c[i]*m[i];
	return r;
}

double eval_pol20_dy(double c[20], double x, double y, double z)
{
	double m[20] = {0, 0, y, 0, x,
		0, z, 0, 2*y, 0,
		x*z, 0, x*2*y, 0, x*x,
		3*y*y, z*z, 0, 2*y*z, 0};
	double r = 0;
	for (int i = 0; i < 20; i++)
		r += c[i]*m[i];
	return r;
}

double eval_pol20_dz(double c[20], double x, double y, double z)
{
	double m[20] = {0, 0, 0, z, 0,
		x, y, 0, 0, 2*z,
		x*y, 0, 0, x*2*z, 0,
		0, y*2*z, x*x, y*y, 3*z*z};
	double r = 0;
	for (int i = 0; i < 20; i++)
		r += c[i]*m[i];
	return r;
}

// evaluate the normalized direct rpc model
static void eval_nrpc(double *result,
		struct rpc *p, double x, double y, double z)
{
	double numx = eval_pol20(p->numx, x, y, z);
	double denx = eval_pol20(p->denx, x, y, z);
	double numy = eval_pol20(p->numy, x, y, z);
	double deny = eval_pol20(p->deny, x, y, z);
	result[0] = numx/denx;
	result[1] = numy/deny;
	//fprintf(stderr, "\t\tnrpc{%p}(%g %g %g)=>(%g %g)\n", p, x, y, z, result[0], result[1]);

}

// evaluate the normalized inverse rpc model
static void eval_nrpci(double *result,
		struct rpc *p, double x, double y, double z)
{
	double numx = eval_pol20(p->inumx, x, y, z);
	double denx = eval_pol20(p->idenx, x, y, z);
	double numy = eval_pol20(p->inumy, x, y, z);
	double deny = eval_pol20(p->ideny, x, y, z);
	result[0] = numx/denx;
	result[1] = numy/deny;
	//fprintf(stderr, "\t\tnrpci{%p}(%g %g %g)=>",p,x,y,z);
	//fprintf(stderr, "(%g %g)\n", result[0], result[1]);
}

// evaluate the direct rpc model
void eval_rpc(double *result,
		struct rpc *p, double x, double y, double z)
{
	double nx = (x - p->offset[0])/p->scale[0];
	double ny = (y - p->offset[1])/p->scale[1];
	double nz = (z - p->offset[2])/p->scale[2];
	double tmp[2];
	eval_nrpc(tmp, p, nx, ny, nz);
	result[0] = tmp[0] * p->iscale[0] + p->ioffset[0];
	result[1] = tmp[1] * p->iscale[1] + p->ioffset[1];
}


// evaluate the inverse rpc model
void eval_rpci(double *result,
		struct rpc *p, double x, double y, double z)
{
	double nx = (x - p->ioffset[0])/p->iscale[0];
	double ny = (y - p->ioffset[1])/p->iscale[1];
	double nz = (z - p->ioffset[2])/p->iscale[2];
	double tmp[2];
	eval_nrpci(tmp, p, nx, ny, nz);
	result[0] = tmp[0] * p->scale[0] + p->offset[0];
	result[1] = tmp[1] * p->scale[1] + p->offset[1];
}

void rpc_projection(double out_pq[2], struct rpc *r,
		double lon, double lat, double h)
{
	eval_rpci(out_pq, r, lon, lat, h);
}

void rpc_localization(double out_latlon[2],
		struct rpc *r, double p, double q, double h)
{
	eval_rpc(out_latlon, r, p, q, h);
}


// evaluate a correspondence between to images given their rpc
void eval_rpc_pair(double xprime[2],
		struct rpc *pa, struct rpc *pb,
		double x, double y, double z)
{
	double tmp[2];
	eval_rpc(tmp, pa, x, y, z);
	eval_rpci(xprime, pb, tmp[0], tmp[1], z);
}

// evaluate the derivative with respect to h of rpc_pair
void eval_rpc_pair_dh(double dph[2],
		struct rpc *pa, struct rpc *pb,
		double x, double y, double h)
{

}

#include "smapa.h"
SMART_PARAMETER(HSTEP,1)
//SMART_PARAMETER(A2TOP,1e-20)
SMART_PARAMETER(LAMAX,1e-6)

// compute the height of a point given its location inside two images
double rpc_height2(struct rpc *rpca, struct rpc *rpcb,
		double xa, double ya, double xb, double yb, double outerr[2])
{
	double x[2] = {xa, ya};
	double y[2] = {xb, yb};

	double h = 0;
	for (int t = 0; t < 1000; t++) {
		double hstep = HSTEP();
		double p[2], q[2];
		eval_rpc_pair(p, rpca, rpcb, x[0], x[1], h);
		eval_rpc_pair(q, rpca, rpcb, x[0], x[1], h + hstep);

		double a[2] = {q[0] - p[0], q[1] - p[1]};
		double b[2] = {y[0] - p[0], y[1] - p[1]};
		double a2 = a[0]*a[0] + a[1]*a[1];

		//if (a2 < A2TOP()) break;

		double lambda = (a[0]*b[0] + a[1]*b[1])/a2;

		// projection of p2 to the line r1-r0
		double z[2] = {p[0] + lambda*a[0], p[1] + lambda*a[1]};

		outerr[0] = z[0] - y[0];
		outerr[1] = z[1] - y[1];

		h += lambda*hstep;

		if (fabs(lambda) < LAMAX())
			break;
	}
	return h;
}

// compute the height of a point given its location inside two images
double rpc_height(struct rpc *rpca, struct rpc *rpcb,
		double xa, double ya, double xb, double yb, double *outerr)
{
	double e[2];
	double r = rpc_height2(rpca, rpcb, xa, ya, xb, yb, e);
	*outerr = hypot(e[0], e[1]);
	return r;
}



static double random_uniform(void)
{
	return rand()/(RAND_MAX+1.0);
}

static int main_trial2(int c, char *v[])
{
	struct rpc p[1];
	nan_rpc(p);
	read_rpc_file_xml(p, "-");
	double lx[][3] = { {1,1,0}, {41500,1,0}, {1,16992,0}, {41500,16992,0}};
	for(int i = 0; i < 4; i++)
	{
		fprintf(stderr, "\n");
		double x = lx[i][0];
		double y = lx[i][1];
		double z = lx[i][2];
		double r[2], rr[2];
		eval_rpc(r, p, x, y, z);
		eval_rpci(rr, p, r[0], r[1], z);
		fprintf(stderr, "(%g %g %g) => ", x, y, z);
		fprintf(stderr, "(%d:%d %d:%d)",
		(int)trunc(r[0]), (int)trunc(60*fabs(r[0]-trunc(r[0]))),
		(int)trunc(r[1]), (int)trunc(60*fabs(r[1]-trunc(r[1])))
				);
		fprintf(stderr, " => (%g %g)\n", rr[0], rr[1]);
	}
	return 0;
}

static int main_trial(int c, char *v[])
{
	struct rpc p[1];
	nan_rpc(p);
	read_rpc_file_xml(p, "-");
	print_rpc(stderr, p, "p");
	for (int i = 0; i < 10; i++) {
		double x = random_uniform();
		double y = random_uniform();
		double z = random_uniform();
		double r[2], rr[2];
		eval_nrpc(r, p, x, y, z);
		eval_nrpci(rr, p, r[0], r[1], z);
		fprintf(stderr, "(%g %g %g) => (%g %g) => (%g %g)\n",
				x, y, z, r[0], r[1], rr[0], rr[1]);
		//fprintf(stderr, "%g\n", hypot(rr[0]-x, rr[1]-y));
	}
	for (int i = 0; i < 10; i++) {
		double x = 10000+4000*random_uniform();
		double y = 10000+4000*random_uniform();
		double z = 1000*random_uniform();
		double r[2], rr[2];
		eval_rpc(r, p, x, y, z);
		eval_rpci(rr, p, r[0], r[1], z);
		fprintf(stderr, "(%g %g %g) => (%g %g) => (%g %g)\n",
				x, y, z, r[0], r[1], rr[0], rr[1]);
		//fprintf(stderr, "%g\n", hypot(rr[0]-x, rr[1]-y));
	}
	return 0;
}

#include "smapa.h"
SMART_PARAMETER(RPCL_NH,21)
SMART_PARAMETER(RPCL_EQUI,0)

static int main_rpcline(int c, char *v[])
{
	if (c != 11) {
		fprintf(stderr, "usage:\n\t"
			"%s a0x a0y b0x b0y rpca rpcb x y h0 hf"
		//        0 1   2   3   4   5    6    7 8 9  10
			"\n", *v);
		return EXIT_FAILURE;
	}
	int offset_a[2] = {atoi(v[1]), atoi(v[2])};
	int offset_b[2] = {atoi(v[3]), atoi(v[4])};
	char *filename_rpca = v[5];
	char *filename_rpcb = v[6];
	double basepoint[2] = {atof(v[7]), atof(v[8])};
	double hrange[2] = {atof(v[9]), atof(v[10])};

	struct rpc rpca[1]; read_rpc_file_xml(rpca, filename_rpca);
	struct rpc rpcb[1]; read_rpc_file_xml(rpcb, filename_rpcb);
	//print_rpc(stderr, rpca, "a");
	//print_rpc(stderr, rpcb, "b");
	int nh = RPCL_NH();

	for (int i = 0; i < nh; i++)
	{
		double ix = basepoint[0];
		double iy = basepoint[1];
		double x = ix + offset_a[0];
		double y = iy + offset_a[1];
		double ni = i/(nh - 1.0);
		if (!RPCL_EQUI()) ni = ni*ni;
		double z = hrange[0] + ni * (hrange[1] - hrange[0]);
		double r[2], rr[2];
		//fprintf(stderr, "(%g %g %g) =>\t", ix, iy, z);
		eval_rpc(r, rpca, x, y, z);
		//fprintf(stderr, "(%.10g %.10g) =>\t", r[0], r[1]);
		eval_rpci(rr, rpcb, r[0], r[1], z);
		double ox = rr[0] - offset_b[0];
		double oy = rr[1] - offset_b[1];
		//fprintf(stderr, "\n\tr[0] = %g\n", r[0]);
		//fprintf(stderr, "\tr[1] = %g\n", r[1]);
		//fprintf(stderr, "\trr[0] = %g\n", rr[0]);
		//fprintf(stderr, "\trr[1] = %g\n", rr[1]);
		//fprintf(stderr, "(%.20g %.20g)\n", ox, oy);
		fprintf(stdout, "%.20g %.20g\n", ox, oy);
	}
	return 0;
}

static int main_rpcpair(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s rpca rpcb x y h\n", *v);
		//                          0 1    2    3 4 5
		return EXIT_FAILURE;
	}
	char *filename_rpca = v[1];
	char *filename_rpcb = v[2];
	double x[3] = {atof(v[3]), atof(v[4]), atof(v[5])}, y[2];

	struct rpc a[1]; read_rpc_file_xml(a, filename_rpca);
	struct rpc b[1]; read_rpc_file_xml(b, filename_rpcb);

	eval_rpc_pair(y, a, b, x[0], x[1], x[2]);

	printf("%.20g %.20g\n", y[0], y[1]);

	return 0;
}


#ifndef DONT_USE_TEST_MAIN
int main(int c, char *v[])
{
	//return main_trial(c, v);
	//return main_trial2(c, v);
	return main_rpcline(c, v);
	//return main_rpcpair(c, v);
}
#endif
