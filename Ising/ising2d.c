/* Program to solve the 2-dimensional Ising model by Metropolis algorithm.
 * Calculates the absolute magnetisation, susceptibility and avgerage energy
 * over a range of temperatures near the phase transition region.
 */
#include <stdio.h>
#include <math.h>
#include "util.h"
#define DIM 32
#define NSPINS (DIM * DIM)

short lat[DIM][DIM];		/* Global lattice variable */
enum lattype { cold, hot };

/* Generate a "hot" (random) or "cold" (all up spins) lattice.
 * Returns the net magnetisation of the resulting lattice */
long genlat(enum lattype type)
{
	int i, j;
	long mag;

	for (mag = i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j) {
			if (type == cold)
				lat[i][j] = 1;
			else
				/* 2 * { 0, 1 } - 1 = { -1, 1 } */
				lat[i][j] = 2 * round(frand(0, 1)) - 1;
			mag += lat[i][j];
		}
	}
	return mag;
}

/* Return the total energy of the lattice, given by the Ising Hamiltonian */
long elat()
{
	int i, j;
	long egy;

	for (egy = i = 0; i < DIM; ++i)
		for (j = 0; j < DIM; ++j)
			egy -= lat[i][j] * (lat[i][(j + 1) % DIM]
					+ lat[(i + 1) % DIM][j]);
	return egy;
}

/* Calculate the energy difference that would be produced by flipping the
 * (i, j)th spin, without actually flipping it */
short eflip(int i, int j)
{
	int n, e, w, s;

	/* Find neighbours, respecting periodic boundary conditions */
	n = (i + 1) % DIM;
	e = (j + 1) % DIM;
	w = (j - 1 + DIM) % DIM;
	s = (i - 1 + DIM) % DIM;
	/* dE = 2s(s1 + s2 + s3 + s4) */
	return 2 * lat[i][j] * (lat[n][j] + lat[i][e] + lat[i][w] + lat[s][j]);
}

/* Calculate the Metropolis probabilities for the given temperature and store
 * them in prob[] for quick access. In 2D Ising model, only two positive values
 * (excluding zero) are possible for delta-E; hence we need to store only two
 * probabilities. */
void metropolis(double temp, double prob[2])
{
	const static short ediffs[] = { 4, 8 };
	int i;

	for (i = 0; i < LEN(ediffs); ++i)
		prob[i] = exp(-ediffs[i] / temp);
}

/* Perform a single Monte-Carlo step (sweep) of the lattice.
 * Precalculated probabilities are to be supplied to p[]. <e> and <m> are the
 * instantaneous energy and magnetisation variables respectivey, to be updated
 * during the sweep (so that we don't have to run a loop to measure them again
 * and again). The average (over the lattice) of quantities will be calculated
 * and returned, if any of the out parameters <egy>, <mag> or <msqr> are given.
 */
void sweep(double p[2], long *e, long *m, double *egy, double *mag, double *msqr)
{
	short ediff;
	int i, j;

	for (i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j) {
			ediff = eflip(i, j);
			/* The Metropolis acceptance criteria */
			/* ediff is converted to array index if > 0 */
			if (ediff <= 0 || frand(0, 1) <= p[ediff / 4 - 1]) {
				lat[i][j] = -lat[i][j];		/* flip */
				/* Update the change in energy & magnetisation */
				*e += ediff;
				*m += 2 * lat[i][j];
			}
			/* Measure, if parameters are given (not NULL) */
			if (egy) *egy += (double) *e / NSPINS;
			if (mag) *mag += fabs((double) *m) / NSPINS;
			if (msqr) *msqr += pow((double) *m / NSPINS, 2);
		}
	}
}

/* Main loop: simulate the 2D Ising model and record the parameters over a range
 * of temperatures. Both adaptive step-sizing (to step through temperatures) and
 * automatic equilibrium detection are used for maximum optimisation. */
int main()
{
	const int steps = 1e5;		/* No. of sweeps */
	int i, k;
	long m, e;
	double temp, egy, mag, msqr, chi, prob[2];
	double de1, de2, olde, oldegy, h;

	/* Initialise: we choose a "hot-start" */
	m = genlat(hot);
	e = olde = elat();
	h = 0.1;	/* Default temperature stepsize */
	k = 2;		/* No. of measurements needed to start adapting h */
	for (temp = 5; temp >= 0; temp -= h) {
		fprintf(stderr, "t = %lf\t", temp);
		metropolis(temp, prob);
		/* Equilibriation sweeps: without measuring magnetisation */
		for (i = 0; i < steps; ++i) {
			egy = 0.0;
			sweep(prob, &e, &m, &egy, NULL, NULL);
			/* Exit early if energy stabilised enough */
			if (eq(egy, olde)) break;
			olde = egy;
		}
		fprintf(stderr, "eq in %d steps\n", i + 1);
		mag = msqr = egy = 0.0;
		/* Measurement sweeps */
		for (i = 0; i < steps; ++i)
			sweep(prob, &e, &m, &egy, &mag, &msqr);

		/* Time averaging */
		mag /= steps * NSPINS;
		msqr /= steps * NSPINS;
		egy /= steps * NSPINS;
		chi = (msqr - mag * mag) / temp;
		printf("%e\t%e\t%e\t%e\n", temp, mag, chi, egy);

		de2 = fabs(egy - oldegy);
		if (--k < 0)
			/* Check variation in energy over adjacent temperature
			 * steps and update stepsize */
			h *= cbrt(de1 / de2);
		oldegy = egy;
		de1 = de2;
	}
	return 0;
}
