/* Program to compute the Lyapunov spectrum of various smooth dynamical systems
 * by the  method of tangent-vectors evolution with periodic Gram-Schmidt ortho-
 * normalization.
 * Language: C	(standard: C99)					Version: 1.0
 * Author: Shivaprasad V					Date: 5 Oct 2019
 * Credentials: PH18C032, M.Sc. Physics '18-'20, IIT Madras, IND
 */
#include <math.h>
#include "util.h"
#include "henon_heiles.h"
/* Lorenz attractor parameters */
#define SIGMA 16.0
#define RHO 45.92
#define BETA 4.0
/* Rössler attractor parameters */
#define A 0.15
#define B 0.20
#define C 10.0

/* Compute the time derivatives of coordinates r[] & store them in rdot[] */
void lorenz(double t, double r[3], double rdot[3])
{
	double x, y, z;
	x = r[0], y = r[1], z = r[2];

	/* Lorenz equations */
	rdot[0] = SIGMA * (y - x);
	rdot[1] = x * (RHO - z) - y;
	rdot[2] = x * y - BETA * z;
}

/* Compute the time derivative of tangents in phase space for Lorenz system.
 * w[] should contain the position vector of the current point, followed by
 * the 3 tangent vectors => 12 components; same order for derivatives in wdot[]
 */
void lortan(double t, double w[12], double wdot[12])
{
	int i;
	double x, y, z, w1, w2, w3;
	x = w[0], y = w[1], z = w[2];	/* store the current x, y, z coord. */

	lorenz(t, w, wdot);		/* evolve the coordinates part */
	for (i = 3; i < 12; i += 3) {	/* loop through each tangent vector */
		w1 = w[i], w2 = w[i + 1], w3 = w[i + 2];
		/* Lorenz system tangent evolution equations */
		wdot[i] = SIGMA * (w2 - w1);
		wdot[i + 1] = (RHO - z) * w1 - w2 - x * w3;
		wdot[i + 2] = x * w2 + y * w1 - BETA * w3;
	}
}

/* Compute the time derivatives of coordinates r[] & store them in rdot[] */
void rossler(double t, double r[3], double rdot[3])
{
	double x, y, z;
	x = r[0], y = r[1], z = r[2];

	/* Rössler attractor equations */
	rdot[0] = -y - z;
	rdot[1] = x + A * y;
	rdot[2] = B + z * (x - C);
}

/* Compute the time derivative of tangents in phase space for Rössler attractor.
 * Identical to lortan() above */
void rosstan(double t, double w[12], double wdot[12])
{
	int i;
	double x, y, z, w1, w2, w3;
	x = w[0], y = w[1], z = w[2];	/* store the current x, y, z coord. */

	rossler(t, w, wdot);         	/* evolve the coordinates part */
	for (i = 3; i < 12; i += 3) {	/* loop through each tangent vector */
		w1 = w[i], w2 = w[i + 1], w3 = w[i + 2];
		/* Rössler attractor tangent evolution equations */
		wdot[i] = -w2 - w3;
		wdot[i + 1] = w1 + A * w2;
		wdot[i + 2] = z * w1 + (x - C) * w3;
	}
}

/* Compute and print the Lyapunov spectrum of a <d> dimensional system whose
 * tangent-evolution is specified by <tanevol>. r[] is the initial position-
 * vector & <b> is the base to which Lyapunov exponents are to be calculated.
 */
void lyaspec(deriv_func tanevol, double r0[], int d, double b)
{
	const int n = 1e7;		/* no. of RK4 steps */
	const int l = (d + 1) * d;	/* no. of components to evolve */
	int i, t;
	double h;
	double w[l], *v[d], a[d], lna[d];

	/* Initialize */
	scale(lna, d, 0.0);		/* set lna to [0 0 ... 0] */
	scale(w, l, 0.0);		/* set w to [0 0 ... 0] */
	copy_vec(r0, w, d);		/* store coordinates at the beg of w */
	for (i = 0; i < d; ++i) {
		v[i] = &w[(i + 1) * d];		/* pointer to i-th tan vector */
		v[i][i] = 1.0;			/* make orthonormal initially */
	}
	h = 0.001;
	t = 0;
	while (t < n) {
		rk4(h * t++, w, l, tanevol, h, w);
		if (t % 100 == 0) {	/* orthonormalise every 100th step */
			ortho(v, d, d, a);
			/* Sum over the natural log of each scaling factor */
			for (i = 0; i < d; ++i)
				lna[i] += log(a[i]);
		}
	}
	/* Divide by total time and convert the exponent to base <b> */
	scale(lna, d, 1.0 / (n * h * log(b)));
	print_vec(lna, d);
}

/* Main program: compute the Lyapunov spectrum (to the base 'e') of 3 systems */
int main()
{
	double egy, r[4];

	/* Lorenz and Rössler attractors */
	r[0] = 1.0; r[1] = 0.0; r[2] = 5.0;
	lyaspec(lortan, r, 3, M_E);
	lyaspec(rosstan, r, 3, M_E);
	/* Hénon-Heiles system: for an invariant trajectory */
	egy = 1 / 6.0;
	r[0] = r[3] = 0.0;
	r[1] = 0.5;
	r[2] = px(egy, r[0], r[1], r[3]);
	lyaspec(hh_tan, r, 4, M_E);
	/* Hénon-Heiles system: for an ergodic trajectory */
	r[1] = 0.3;
	r[2] = px(egy, r[0], r[1], r[3]);
	lyaspec(hh_tan, r, 4, M_E);
	return 0;
}

/* Extra functions: */
/* Plot the Lorenz attractor for <n> RK4 steps */
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

/* Plot the Rössler attractor for <n> RK4 steps */
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
