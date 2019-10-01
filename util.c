/* Utilities for genral purpose programming */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util.h"

/* Dynamically allocate an n element array of doubles */
double *vector(int n)
{
	double *v;

	v = (double *) malloc((size_t) n * sizeof(double));
	return v;
}

/* Compute inner-product of two vectors: u.v */
double dot(double *u, double *v, int n)
{
	double s;

	for (s = 0.0; n > 0; --n)
		s += *u++ * *v++;
	return s;
}

/* Multiply vector by a scalar k */
void scale(double *v, int n, double k)
{
	for (; n > 0; --n)
		*v++ *= k;
}

/* Return the projection of v along u as a dynamically allocated vector *
 * proj(v, u) = v.u * u / u.u */
double *proj(double *v, double *u, int n)
{
	double *t;

	t = vector(n);
	copy_vec(u, t, n);
	scale(t, n, dot(u, v, n) / dot(u, u, n));
	return t;
}

/* Subtract vector u from v: v - u */
void sub(double *v, double *u, int n)
{
	for (; n > 0; --n)
		*v++ -= *u++;
}

/* Normalise given vector v */
void norm(double *v, int n)
{
	double a;

	a = sqrt(dot(v, v, n));
	for (; n > 0; --n)
		*v++ /= a;
}

/* Print the elements of an array of size n */
void print_vec(double *v, int n)
{
	for (; n > 0; --n)
		printf("%.7lf\t", *v++);
	putchar('\n');
}

/* Copy array s to t, both having size n */
void copy_vec(double *s, double *t, int n)
{
	for (; n > 0; --n)
		*t++ = *s++;
}

/* Calculate the n dimensional distance between two vectors */
double dist(double *v1, double *v2, int n)
{
	double d;

	for (d = 0.0; n > 0; --n)
		d += pow((*v1++ - *v2++), 2);
	return sqrt(d);
}

/* Return a random real number between lb and ub (both inclusive).
 * Employs the Bays-Durham shuffling algorithm as given in Num. Recep. */
double frand(double lb, double ub)
{
	int i, r;
	static int p = -1;		/* next rand no. chooser */
	static int deck[DECK_SIZE];	/* shuffling table */

	if (p < 0) {
		srand(time(NULL));
		/* Load the table after 8 warmups */
		for (i = -8; i < DECK_SIZE; ++i) {
			r = rand();
			if (i >= 0) deck[i] = r;
		}
		p = rand();
	}
	i = (int) (p / (RAND_MAX + 1.0) * DECK_SIZE);	/* i = 0...31 */
	p = deck[i];		/* next chooser */
	deck[i] = rand();	/* refill vacancy */
	return (ub - lb) * p / RAND_MAX + lb;	/* scale & return */
}

/* Compare two real numbers; tolerance can be defined as EPS */
enum bool eq(double x, double y)
{
	if (fabs(x - y) < EPS)
		return true;
	else
		return false;
}

/* n dimensional RK4 algorithm for calculating a single step.
 * This implementation is loosely based on the one given in Numerical Recipes.
 * As we cannot multiply arrays with numbers in C (without using a loop), an
 * approach different from the traditional way (k1, k2, ...) has been followed.
 */
void rk4(double x, double y[], int n, deriv_func f, double h, double yout[])
{
	int i;
	double hh, xh, *dy, *dys, *yt;

	yt = vector(n);
	dy = vector(n);		/* derivatives of y1, y2, ... */
	dys = vector(n);	/* sum of derivatives */
	hh = 0.5 * h;
	xh = x + hh;

	f(x, y, dy);				/* first step */
	for (i = 0; i < n; ++i) {
		yt[i] = y[i] + hh * dy[i];	/* y + k1 / 2 */
		dys[i] = dy[i];
	}
	f(xh, yt, dy);				/* second step */
	for (i = 0; i < n; ++i) {
		yt[i] = y[i] + hh * dy[i];	/* y + k2 / 2 */
		dys[i] += 2 * dy[i];
	}
	f(xh, yt, dy);				/* third step */
	for (i = 0; i < n; ++i) {
		yt[i] = y[i] + h * dy[i];	/* y + k3 */
		dys[i] += 2 * dy[i];
	}
	f(x + h, yt, dy);			/* fourth step */
	for (i = 0; i < n; ++i)
		yout[i] = y[i] + h / 6.0 * (dys[i] + dy[i]);	/* final y */

	free(yt);
	free(dy);
	free(dys);
}

/* Orthonormalise a set of m n-dimensional vectors by Gram–Schmidt procedure */
void ortho(double **u, int m, int n)
{
	int i, k;
	double *p;

	for (k = 0; k < m; ++k) {
		/* Subtract components along all u[i < k] from u[k] */
		for (i = 0; i < k; ++i) {
			p = proj(u[k], u[i], n);
			sub(u[k], p, n);
			free(p);
		}
		norm(u[k], n);	/* normalise */
	}
}
