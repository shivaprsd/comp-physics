#include <math.h>
#include "util.h"
#define DIM 3
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

void lortang(double t, double w[12], double wdot[12])
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

void rosstang(double t, double w[12], double wdot[12])
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

int main()
{
	const int n = 1e7;
	const int l = (DIM + 1) * DIM;
	int i, t;
	double h;
	double w[l], *v[DIM], a[DIM], lna[DIM];

	scale(lna, DIM, 0.0);
	scale(w, l, 0.0);
	w[0] = 1.0; w[1] = 0.0; w[2] = 5.0;
	for (i = 0; i < DIM; ++i) {
		v[i] = &w[(i + 1) * DIM];
		v[i][i] = 1.0;
	}
	h = 0.001;
	t = 0;
	while (t < n) {
		rk4(h * t++, w, l, lortang, h, w);
		if (t % 100 == 0) {
			ortho(v, DIM, DIM, a);
			for (i = 0; i < DIM; ++i)
				lna[i] += log(a[i]);
		}
	}
	scale(lna, DIM, 1.0 / (n * h));
	print_vec(lna, DIM);
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
