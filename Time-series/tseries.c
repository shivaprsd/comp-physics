/* Program to make short-term predictions of the future behaviour of a chaotic
 * time series from known data by local-approximation method outlined by Farmer
 * and Sidorowich in Phys. Rev. Lett. #59 (1987) 845-48.
 * Language: C	(standard: ANSI C)				Version: 1.0
 * Author: Shivaprasad V					Date: 17 Nov 2019
 * Credentials: PH18C032, M.Sc. Physics '18-'20, IIT Madras, IND
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>		/* for DBL_MAX */
#include "util.h"

#define NPRED 4096
#define NMAX (12 * NPRED)	/* Reasonable sample size for given NPRED */
#define TMAX 15
/* Parameters for Logistic and Henon maps */
#define MU 4.0
#define A 1.4
#define B 0.3
enum dynsys { LOG_MAP, HEN_MAP, REAL_WLD };

/* Generate (or read) <n> terms (including initial points) of a time series of
 * given <type> and store them in <s> */
void genseries(double *s, int n, enum dynsys type)
{
	int i;
	double x, y;

	switch (type) {
		case LOG_MAP:
			for (i = 1; i < n; ++i) {	/* one initial point */
				x = s[i - 1];
				s[i] = MU * x * (1 - x);	/* Log. map */
			}
			break;
		case HEN_MAP:
			/* We need only a one-variable series; y is discarded */
			for (i = 2; i < n; ++i) {	/* two initial points */
				x = s[i - 1];
				y = B * s[i - 2];
				s[i] = 1 - A * x * x + y;	/* Henon map */
			}
			break;
		case REAL_WLD:
			/* Feed any other (real-world) series through stdin */
			for (i = 0; i < n; ++i)
				scanf("%lf\n", s + i);
			break;
	}
}

/* Do a linear least-square fit of <n> pairs of elements in given domain and
 * range by QR decompostion. <d> is the dimension of each element of <dom> which
 * can be vectors, while <rng> should be a one-dimensional array of <n> scalars.
 * Returns the fit-parameters b[i] as a dynamically allocated array. */
double *linfit(double **dom, int n, int d, double *rng)
{
	int i, j;
	double **a, **q, **r, *b;

	/* We need d + 1 locations for the d components plus one constant term;
	 * So increment the local variable to avoid writing d + 1 everywhere */
	a = matrix(++d, n);
	q = matrix(d, n);
	r = matrix(d, d);	/* assuming n > d */
	set(a[0], n, 1.0);	/* column of coeff. of the constant term = 1 */
	for (j = 1; j < d; ++j) {
		for (i = 0; i < n; ++i)
			/* store the components in column-major order */
			a[j][i] = dom[i][j - 1];
	}
	qrd(a, q, r, n, d);	/* QR decompose A */
	b = vector(d);
	/* Need to solve Rx = Q'b; as Q is stored column-wise, Q'b = dot(Q,b) */
	for (i = 0; i < d; ++i)
		b[i] = dot(q[i], rng, n);
	bsub(r, d, b);		/* solve in-place to obtain fit parameters */

	/* Cleanup allocated memory and return b */
	for (j = 0; j < d; ++j) {
		free(a[j]);
		free(q[j]);
		free(r[j]);
	}
	free(a); free(q); free(r);
	return b;
}

/* Predict next <tmax> terms of time series <data> having <n> elements and store
 * them in <pred>. <dim> is the embedding dimension to use & <ord> is the number
 * of nearest neighbours to consider (= 1 + order of approximation). */
void predict(double *data, int n, int dim, int ord, int tmax, double *pred)
{
	int i, j, t, *nbr;
	double d, *dmin, *x, *xt, *fp, **p;

	/* Need more points than parameters for convergence of linear fit */
	if (ord != 1 && ord <= dim) {
		fprintf(stderr, "Error: need to consider %d or more neighbours "
				"for a %d-dimensional system!\n", dim + 1, dim);
		return;
	}
	nbr = (int *) malloc((size_t) ord * sizeof(int));
	dmin = vector(ord);
	set(dmin, ord, DBL_MAX);	/* for initial comparisons (below) */
	x = data + n - dim;		/* beginning of the reference vector */
	/* Last <tmax> points are not available for 'learning', as x[t' + tmax]
	 * goes out of range; last <dim> points are components of the reference
	 * vector and can't be considered as well; mark off the bigger chunk */
	n -= (dim > tmax) ? dim : tmax;

	/* Find the starting indices of the nearest-neighbour vectors */
	for (i = 0; i < n; ++i) {
		d = dist(data + i, x, dim);
		for (j = 0; j < ord; ++j) {		/* insertion sort */
			if (d < dmin[j]) {
				fins(d, dmin, ord, j);
				ins(i, nbr, ord, j);
				break;
			}
		}
	}
	/* 0th order approximation */
	if (ord == 1) {
		for (t = 0; t < tmax; ++t)
			/* t' = neighbour's starting index + vector length - 1 */
			pred[t] = data[nbr[0] + dim + t];
	/* 1st order approximation */
	} else {
		p = matrix(ord, 0);
		xt = dmin;	/* just reusing previously allocated vector */
		for (j = 0; j < ord; ++j)
			p[j] = data + nbr[j];
		for (t = 0; t < tmax; ++t) {
			for (j = 0; j < ord; ++j)
				xt[j] = *(p[j] + dim + t);
			/* Do a linear fit of x[t'] and x[t' + T]; fp[0] is the
			 * constant coefficient & the rest are coeff. of x[t] */
			fp = linfit(p, ord, dim, xt);
			pred[t] = fp[0] + dot(x, fp + 1, dim);
		}
		free(p);
		free(fp);
	}
	free(nbr);
	free(dmin);
}

/* Main loop: to calculate the normalised error in the predictions for the range
 * of prediction time [1, T]; demonstrating for Logistic map, 1st order approx.
 */
int main()
{
	int i, j, t;
	const int d = 1, k = 3;
	const int ndat = NMAX - (NPRED + TMAX - 1); 	/* no. of data points */
	double x[NMAX], px[TMAX], err[TMAX];

	/* x[0] = 0.1; x[1] = 1.086; can be taken for Henon map */
	x[0] = 3.0 / 8;
	genseries(x, NMAX, LOG_MAP);
	set(err, TMAX, 0.0);
	for (i = 0; i < NPRED; ++i) {			/* start predicting */
		t = ndat + i;
		predict(x, t, d, k, TMAX, px);
		for (j = 0; j < TMAX; ++j)
			err[j] += pow(px[j] - x[t + j], 2);
	}
	for (j = 0; j < TMAX; ++j)			/* normalisation */
		err[j] = sqrt(err[j] / NPRED) / rmsd(x + ndat + j, NPRED);
	printf( "# Time series prediction - normalised errors (for T = 1 to %d)"
		"\n# System: logistic map; embedding dimension: %d; order of "
		"approximation: %d\n", TMAX, d, k );
	print_vec(err, TMAX, "%.7e\n");			/* final output */
	return 0;
}

/* End of program */
