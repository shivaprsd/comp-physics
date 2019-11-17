#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>	/* for DBL_MAX */
#include "util.h"

#define NPRED 4096
#define NMAX (12 * NPRED)
#define TMAX 15
#define MU 4.0
#define A 1.4
#define B 0.3
enum dynsys { LOG_MAP, HEN_MAP, REAL_WLD };

void genseries(double *s, int n, enum dynsys type)
{
	int i;
	double x, y;

	switch (type) {
		case LOG_MAP:
			for (i = 1; i < n; ++i) {
				x = s[i - 1];
				s[i] = MU * x * (1 - x);
			}
			break;
		case HEN_MAP:
			for (i = 2; i < n; ++i) {
				x = s[i - 1];
				y = B * s[i - 2];
				s[i] = 1 - A * x * x + y;
			}
			break;
		case REAL_WLD:
			for (i = 0; i < n; ++i)
				scanf("%lf\n", s + i);
			break;
	}
}

double *linfit(double **dom, int n, int d, double *rng)
{
	int i, j;
	double **a, **q, **r, *b;

	a = matrix(++d, n);
	q = matrix(d, n);
	r = matrix(d, d);
	set(a[0], n, 1.0);
	for (j = 1; j < d; ++j) {
		for (i = 0; i < n; ++i)
			a[j][i] = dom[i][j - 1];
	}
	qrd(a, q, r, n, d);
	b = vector(d);
	for (i = 0; i < d; ++i)
		b[i] = dot(q[i], rng, n);
	bsub(r, d, b);

	for (j = 0; j < d; ++j) {
		free(a[j]);
		free(q[j]);
		free(r[j]);
	}
	free(a); free(q); free(r);
	return b;
}

void predict(double *data, int n, int dim, int ord, int tmax, double *pred)
{
	int i, j, t, *nbr;
	double d, *dmin, *x, *xt, *fp, **p;

	if (ord != 1 && ord <= dim) {
		fprintf(stderr, "Error: need to consider %d or more neighbours "
				"for a %d-dimensional system!\n", dim + 1, dim);
		return;
	}
	nbr = (int *) malloc((size_t) ord * sizeof(int));
	dmin = vector(ord);
	set(dmin, ord, DBL_MAX);
	x = data + n - dim;
	n -= (dim > tmax) ? dim : tmax;

	for (i = 0; i < n; ++i) {
		d = dist(data + i, x, dim);
		for (j = 0; j < ord; ++j) {
			if (d < dmin[j]) {
				fins(d, dmin, ord, j);
				ins(i, nbr, ord, j);
				break;
			}
		}
	}
	if (ord == 1) {
		for (t = 1; t <= tmax; ++t)
			pred[t - 1] = data[nbr[0] + dim - 1 + t];
	} else {
		p = matrix(ord, 0);
		xt = dmin;
		for (j = 0; j < ord; ++j)
			p[j] = data + nbr[j];
		for (t = 1; t <= tmax; ++t) {
			for (j = 0; j < ord; ++j)
				xt[j] = *(p[j] + dim - 1 + t);
			fp = linfit(p, ord, dim, xt);
			pred[t - 1] = fp[0] + dot(x, fp + 1, dim);
		}
		free(p);
		free(fp);
	}
	free(nbr);
	free(dmin);
}

int main()
{
	int i, j, t;
	const int d = 1, k = 3;
	const int ndat = NMAX - (NPRED + TMAX - 1);
	double x[NMAX], px[TMAX], err[TMAX];

	x[0] = 3.0 / 8;
	//x[1] = 1.086;
	genseries(x, NMAX, LOG_MAP);
	set(err, TMAX, 0.0);
	for (i = 0; i < NPRED; ++i) {
		t = ndat + i;
		predict(x, t, d, k, TMAX, px);
		for (j = 0; j < TMAX; ++j)
			err[j] += pow(px[j] - x[t + j], 2);
	}
	for (j = 0; j < TMAX; ++j)
		err[j] = sqrt(err[j] / NPRED) / rmsd(x + ndat + j, NPRED);
	printf( "# Time series prediction - normalised errors (for T = 1 to %d)"
		"\n# System: logistic map; embedding dimension: %d; order of "
		"approximation: %d\n", TMAX, d, k );
	print_vec(err, TMAX, "%.7e\n");
	return 0;
}
