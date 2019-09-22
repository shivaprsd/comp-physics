#include <stdio.h>
#include <math.h>
#include "util.h"
#define DIM 32
#define NSPINS (DIM * DIM)

short lat[DIM][DIM];
enum lattype { cold, hot };

long genlat(enum lattype type)
{
	int i, j;
	long mag;

	for (mag = i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j) {
			if (type == cold)
				lat[i][j] = 1;
			else
				lat[i][j] = 2 * round(frand(0, 1)) - 1;
			mag += lat[i][j];
		}
	}
	return mag;
}

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

short eflip(int i, int j)
{
	int n, e, w, s;

	n = (i + 1) % DIM;
	e = (j + 1) % DIM;
	w = (j - 1 + DIM) % DIM;
	s = (i - 1 + DIM) % DIM;
	return 2 * lat[i][j] * (lat[n][j] + lat[i][e] + lat[i][w] + lat[s][j]);
}

void metropolis(double temp, double prob[2])
{
	const static short ediffs[] = { 4, 8 };
	int i;

	for (i = 0; i < LEN(ediffs); ++i)
		prob[i] = exp(-ediffs[i] / temp);
}

void sweep(double p[2], long *e, long *m, double *egy, double *mag, double *msqr)
{
	short ediff;
	int i, j;

	for (i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j) {
			ediff = eflip(i, j);
			if (ediff <= 0 || frand(0, 1) <= p[ediff / 4 - 1]) {
				lat[i][j] = -lat[i][j];
				*e += ediff;
				*m += 2 * lat[i][j];
			}
			if (egy) *egy += (double) *e / NSPINS;
			if (mag) *mag += fabs((double) *m) / NSPINS;
			if (msqr) *msqr += pow((double) *m / NSPINS, 2);
		}
	}
}

int main()
{
	const int steps = 1e5;
	int i, k;
	long m, e;
	double temp, egy, mag, msqr, chi, prob[2];
	double de1, de2, olde, oldegy, h;

	m = genlat(hot);
	e = olde = elat();
	h = 0.1;
	k = 2;
	for (temp = 5; temp >= 0; temp -= h) {
		fprintf(stderr, "t = %lf\t", temp);
		metropolis(temp, prob);
		for (i = 0; i < steps; ++i) {
			egy = 0.0;
			sweep(prob, &e, &m, &egy, NULL, NULL);
			if (eq(egy, olde)) break;
			olde = egy;
		}
		fprintf(stderr, "eq in %d steps\n", i + 1);
		mag = msqr = egy = 0.0;
		for (i = 0; i < steps; ++i)
			sweep(prob, &e, &m, &egy, &mag, &msqr);

		mag /= steps * NSPINS;
		msqr /= steps * NSPINS;
		egy /= steps * NSPINS;
		chi = (msqr - mag * mag) / temp;
		printf("%e\t%e\t%e\t%e\n", temp, mag, chi, egy);

		de2 = fabs(egy - oldegy);
		if (--k < 0)
			h *= cbrt(de1 / de2);
		oldegy = egy;
		de1 = de2;
	}
	return 0;
}

/*
void printlat()
{
	int i, j;

	for (i = 0; i < DIM; ++i) {
		putchar('\t');
		for (j = 0; j < DIM; ++j) {
			if (lat[i][j] == -1)
				printf("%lf\t%lf\n", 0.001 * i, 0.001 * j);
		}
		putchar('\n');
	}
}

short flip(double prob[3])
{
	int i, j;
	short ediff;

	i = (int) frand(0, DIM - EPS);
	j = (int) frand(0, DIM - EPS);
	ediff = eflip(i, j);
	if (ediff < 0 || frand(0, 1) <= prob[ediff / 4]) {
		lat[i][j] = -lat[i][j];
		return 2 * lat[i][j];
	}
	return 0;
}

for (i = 0; i < steps * NSPINS; ++i)
	flip(prob);
mag = msqr = 0.0;
for (i = 0; i < steps * NSPINS; ++i) {
	m += flip(prob);
	mag += (double) m / NSPINS;
	msqr += (double) m / NSPINS * (double) m / NSPINS;
}

genlat(cold);
metropolis(2.0, prob);
for (i = 0; i < 200; ++i) {
	sweep(2.0, prob, NULL, NULL);
	printf("%d\t%lf\n", i, (double) elat() / NSPINS);
}
*/
