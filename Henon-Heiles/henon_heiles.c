/* Hénon-Heiles related functions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util.h"

/* Return the value of the Hénon-Heiles potential at (x, y) */
double hh(double x, double y)
{
	return (x * x + y * y) / 2.0 + x * x * y - pow(y, 3) / 3.0;
}

/* Compute the time derivative of v and store it in vdot */
void hh_evol(double t, double v[4], double vdot[4])
{
	double x, y;
	x = v[0], y = v[1];

	/* Hénon-Heiles differential equations of motion */
	vdot[0] = v[2];
	vdot[1] = v[3];
	vdot[2] = -x - 2 * x * y;
	vdot[3] = -y - x * x + y * y;
}

void hh_tan(double t, double w[20], double wdot[20])
{
	int i;
	double x, y, w1, w2;
	x = w[0], y = w[1];

	hh_evol(t, w, wdot);
	for (i = 4; i < 20; i += 4) {
		w1 = w[i], w2 = w[i + 1];
		wdot[i] = w[i + 2];
		wdot[i + 1] = w[i + 3];
		wdot[i + 2] = -(2 * y + 1) * w1 - 2 * x * w2;
		wdot[i + 3] = -2 * x * w1 + (2 * y - 1) * w2;
	}
}

/* Return x-coordinate of the equipotential line u at y
 * Returns -1 if the point is invalid */
double hh_equi(double y, double u)
{
	double xsqr;

	xsqr = (u - y * y / 2 + pow(y, 3) / 3) / (y + 0.5);
	if (xsqr > 0.0 || eq(xsqr, 0.0))
		return sqrt(xsqr);
	else
		return -1;
}

/* Return x-momentum for a given energy, (x, y) and py
 * Returns -1 if the point is invalid */
double px(double e, double x, double y, double py)
{
	double pxsqr;

	pxsqr = 2 * (e - hh(x,y)) - py * py;
	if (pxsqr > 0.0 || eq(pxsqr, 0.0))
		return sqrt(pxsqr);
	else
		return -1;
}

/* Compute the energy from the Hamiltonian, given coordinates x */
double en(double *x, int n)
{
	double ke;

	ke = 0.5 * (pow(x[2], 2) + pow(x[3], 2));
	return ke + hh(x[0], x[1]);
}

/* Compute the next iteration of the area preserving map studied by
 * Hénon and Heiles in their paper, for a given 'a' */
void hh_map(double xin[2], double xout[2], double a)
{
	double xt;

	xt = xin[0] + a * (xin[1] - pow(xin[1], 3));
	xout[1] = xin[1] - a * (xt - pow(xt, 3));
	xout[0] = xt;
}

void hh_map_tan(double win[6], double wout[6], double a)
{
	int i;
	double xn, y, w1, w2;
	y = win[1];

	hh_map(win, wout, a);
	xn = wout[0];
	for (i = 2; i < 6; i += 2) {
		w1 = win[i], w2 = win[i + 1];
		wout[i] = w1 + a * (1 - 3 * y * y) * w2;
		wout[i + 1] = w2 - a * (1 - 3 * xn * xn) * wout[i];
	}
}

/* Print the coordinates of the equipotential lines for a number of u values
 * Prints only the positive half (x > 0) of the curves; mirror to get full */
void plot_equi()
{
	double x, y, u[] = { 1 / 6.0, 1 / 8.0, 1 / 12.0, 1 / 24.0, 1 / 100.0 };
	double xf, yf;
	int i;

	for (i = 0; i < LEN(u); ++i) {
		for (y = 1; y >= -0.5; y -= 0.0001) {
			x = hh_equi(y, u[i]);
			if (x >= 0) {
				printf("%lf\t%lf\n", x, y);
				xf = x;
				yf = y;
			}
		}
		/* Print the mirror of the last point, to get continuation while
		 * plotting */
		printf("%lf\t%lf\n", -xf, yf);
		printf("\n\n");
	}
}

/* Iterate hh_map() n times starting from the initial point x, for given 'a' */
void plot_map(double *x, double a, int n)
{
	int i;

	for (i = 0; i < n; ++i) {
		print_vec(x, 2);
		hh_map(x, x, a);
	}
}

/* Print the first n transformations of the initial point x in the x = 0, px > 0
 * section for a given energy */
void phase_space(double *x, double energy, int n)
{
	int pcount;
	long t;
	double xnew[4], h;

	pcount = t = 0;
	h = 0.001;
	while (pcount < n) {
		rk4(h * t++, x, 4, hh_evol, h, xnew);
		/* Check if we crossed the x = 0 plane with positive momentum */
		if (xnew[0] * x[0] < 0 && x[2] > 0.0) {
			printf("%lf\t%.12lf\t", h * t, en(x, 4));
			print_vec(x, 4);
			++pcount;
		}
		copy_vec(xnew, x, 4);
	}
}

/* End of file */
