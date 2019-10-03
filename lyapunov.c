#include <math.h>
#include "util.h"
/* Lorenz attractor parameters */
#define SIGMA 10.0
#define RHO 28.0
#define BETA (8.0 / 3)
/* Rössler attractor parameters */
#define A 0.15
#define B 0.2
#define C 10.0

void lorenz(double t, double r[3], double rdot[3])
{
	double x, y, z;
	x = r[0], y = r[1], z = r[2];

	rdot[0] = SIGMA * (y - x);
	rdot[1] = x * (RHO - z) - y;
	rdot[2] = x * y - BETA * z;
}

void lortan(double t, double w[12], double wdot[12])
{
	int i;
	double x, y, z, w1, w2, w3;
	x = w[0], y = w[1], z = w[2];

	lorenz(t, w, wdot);
	for (i = 3; i < 12; i += 3) {
		w1 = w[i], w2 = w[i + 1], w3 = w[i + 2];
		wdot[i] = SIGMA * (w2 - w1);
		wdot[i + 1] = (RHO - z) * w1 - w2 - x * w3;
		wdot[i + 2] = x * w2 + y * w1 - BETA * w3;
	}
}

void rossler(double t, double r[3], double rdot[3])
{
	double x, y, z;
	x = r[0], y = r[1], z = r[2];

	rdot[0] = -y - z;
	rdot[1] = x + A * y;
	rdot[2] = B + z * (x - C);
}

void rosstan(double t, double w[12], double wdot[12])
{
	int i;
	double x, y, z, w1, w2, w3;
	x = w[0], y = w[1], z = w[2];

	rossler(t, w, wdot);
	for (i = 3; i < 12; i += 3) {
		w1 = w[i], w2 = w[i + 1], w3 = w[i + 2];
		wdot[i] = -w2 - w3;
		wdot[i + 1] = w1 + A * w2;
		wdot[i + 2] = z * w1 + (x - C) * w3;
	}
}

void lyaspec(deriv_func tanevol, double r0[], int d, double b)
{
	const int n = 1e7;
	const int l = (d + 1) * d;
	int i, t;
	double h;
	double w[l], *v[d], a[d], lna[d];

	scale(lna, d, 0.0);
	scale(w, l, 0.0);
	copy_vec(r0, w, d);
	for (i = 0; i < d; ++i) {
		v[i] = &w[(i + 1) * d];
		v[i][i] = 1.0;
	}
	h = 0.001;
	t = 0;
	while (t < n) {
		rk4(h * t++, w, l, tanevol, h, w);
		if (t % 100 == 0) {
			ortho(v, d, d, a);
			for (i = 0; i < d; ++i)
				lna[i] += log(a[i]);
		}
	}
	scale(lna, d, 1.0 / (n * h * log(b)));
	print_vec(lna, d);
}

int main()
{
	double r[3];

	r[0] = 1.0; r[1] = 0.0; r[2] = 5.0;
	lyaspec(lortan, r, 3, M_E);
	lyaspec(rosstan, r, 3, M_E);
	return 0;
}

void plotlor(int n)
{
	double h, r[3];
	int t;

	r[0] = 1.0; r[1] = 0.0; r[2] = 5.0;
	h = 0.001;
	t = 0;
	while (t < n) {
		rk4(h * t++, r, 3, lorenz, h, r);
		print_vec(r, 3);
	}
}

void plotross(int n)
{
	double h, r[3];
	int t;

	scale(r, 3, 0.0);
	h = 0.001;
	t = 0;
	while (t < n) {
		rk4(h * t++, r, 3, rossler, h, r);
		print_vec(r, 3);
	}
}
