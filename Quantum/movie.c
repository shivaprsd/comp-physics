#include <stdio.h>
#include <math.h>
#include <complex.h>
#define L 1
#define LAMBDA 1
#define V0 (-2 * SQR(50 * M_PI))
#define VW 0.032
#define SQR(x) (pow(x, 2))

int main()
{
	const double sigma_sqr = SQR(L / 20.0);
	const double x0 = L / 4.0;
	const double k0 = 100 * M_PI;
	const int jmax = 20 * k0 / M_PI;
	const int nmax = 1520;
	const double dx = 1.0 / jmax;
	//const double dt = 2 * dx * dx / LAMBDA;
	double x;
	int j, n;

	double v[jmax + 1];
	complex double psi[jmax + 1], omega[jmax], auxe[jmax], auxf[jmax];

	for (j = 0; j <= jmax; ++j) {
		x = j * dx;
		v[j] = (fabs(x - L / 2.0) < VW) ? V0 : 0.0;
		psi[j] = cexp(I * k0 * x) * exp(-SQR(x - x0) / (2 * sigma_sqr));
	}
	//psi[0] = psi[jmax - 1] = 0.0;
	for (j = 1; j < jmax; ++j) {
		auxe[j] = 2 + SQR(dx) * v[j] - I * LAMBDA;
	       	if (j > 1) auxe[j] -= 1 / auxe[j - 1];
	}

	for (n = 0; n <= nmax; ++n) {
		for (j = 1; j < jmax - 1; ++j) {
			omega[j] = (I * LAMBDA + SQR(dx) * v[j] + 2) * psi[j] \
					- psi[j + 1] - psi[j - 1];
			auxf[j] = omega[j];
			if (j > 1) auxf[j] += auxf[j - 1] / auxe[j - 1];
		}
		for (j = jmax - 1; j > 0; --j)
			psi[j] = (psi[j + 1] - auxf[j]) / auxe[j];
		for (j = 0; j <= jmax; ++j)
			printf("%e\t%e\n", SQR(cabs(psi[j])), v[j] / fabs(V0));
		printf("\n\n");
	}
	return 0;
}
