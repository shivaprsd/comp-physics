/* Program to simulate the scattering of a particle having a Gaussian wavepacket
 * by a finite square well / barrier using the method given by Goldberg, Schey &
 * Schwartz in Am. J. Phy. (vol. 35, #3, 1967).
 *
 * Language: C	(standard: C99)					Version: 1.0
 * Author: Shivaprasad V					Date: 1 Nov 2019
 * Credentials: PH18C032, M.Sc. Physics '18-'20, IIT Madras, IND
 */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#define L 1		/* Length of the box */
#define LAMBDA 1	/* 2 (dx)^2 / dt */
#define VW 0.032	/* Half-width of the well / barrier */
#define SQR(x) (pow(x, 2))

/* Main loop: demonstrating for a well with avg. particle energy = well depth */
int main()
{
	/* Define parameters */
	const double sigma_sqr = SQR(L / 20.0);		/* spread of psi */
	const double x0 = L / 4.0;			/* initial peak pos. */
	const double k0 = 50 * sqrt(2) * M_PI;		/* initial momentum */
	const double V0 = -SQR(k0);			/* depth of the well */
	const int jmax = 20 * k0 / M_PI;		/* total mesh points */
	const int nmax = 1000;			/* total simulation timesteps */
	const double dx = 1.0 / jmax;		/* mesh width */

	/* Define variables and mesh arrays */
	int j, n;
	double x;
	double v[jmax + 1];			/* potential values (real) */
	complex double psi[jmax + 1], omega[jmax], auxe[jmax], auxf[jmax];

	/* Store the potential and initial wavefunction at each mesh point */
	for (j = 0; j <= jmax; ++j) {
		x = j * dx;
		psi[j] = cexp(I * k0 * x) * exp(-SQR(x - x0) / (2 * sigma_sqr));
		/* Potential is centered in the box and non-zero only within
		 * a distance of VW away from the center */
		v[j] = (fabs(x - L / 2.0) < VW) ? V0 : 0.0;
		/* Print the (scaled) values of the potential */
		printf("%e\n", v[j] / fabs(V0));
	}
	psi[0] = psi[jmax] = 0.0;		/* boundary conditions */

	/* Store the values of the time-independent auxilary function 'e' */
	for (j = 1; j < jmax; ++j) {
		auxe[j] = 2 + SQR(dx) * v[j] - I * LAMBDA;
		if (j > 1) auxe[j] -= 1 / auxe[j - 1];
	}

	/* Run the simulation for nmax timesteps */
	for (n = 0; n <= nmax; ++n) {
		/* Compute the time-dependent auxilary functions omega & f */
		for (j = 1; j < jmax - 1; ++j) {
			omega[j] = (I * LAMBDA + SQR(dx) * v[j] + 2) * psi[j] \
					- psi[j + 1] - psi[j - 1];
			auxf[j] = omega[j];
			if (j > 1) auxf[j] += auxf[j - 1] / auxe[j - 1];
		}
		/* Time evolve psi by recursing from right to left */
		for (j = jmax - 1; j > 0; --j)
			psi[j] = (psi[j + 1] - auxf[j]) / auxe[j];

		/* Print the probability density |psi|^2 at all mesh points */
		printf("\n\n");
		for (j = 0; j <= jmax; ++j)
			printf("%e\n", SQR(cabs(psi[j])));
	}
	return 0;
}
