#include <stdio.h>
#include <math.h>
#include "util.h"
#define DIM 32
#define NSPINS (DIM * DIM)

struct Site
{
	unsigned int x;
	unsigned int y;
};
char lat[DIM][DIM];

void genlat()
{
	int i, j;

	for (i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j)
			lat[i][j] = 2 * round(frand(0, 1)) - 1;
	}
}

void printlat()
{
	int i, j;

	for (i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j) {
			if (lat[i][j] == -1)
				printf("%lf\t%lf\n", 0.001 * i, 0.001 * j);
		}
	}
}

int elat()
{
	int egy, i, j;

	for (egy = i = 0; i < DIM; ++i)
		for (j = 0; j < DIM; ++j)
			egy -= lat[i][j] * (lat[i][(j + 1) % DIM]
					+ lat[(i + 1) % DIM][j]);
	return egy;
}

long measure()
{
	int i, j;
	long mag;

	mag = 0;
	for (i = 0; i < DIM; ++i) {
		for (j = 0; j < DIM; ++j)
			mag += lat[i][j];
	}
	return mag;
}

void wolff(struct Site *s, double p, long *m)
{
	struct Site nbr[4] = {
		{ (s->x + 1) % DIM, s->y },
		{ s->x, (s->y + 1) % DIM },
		{ (s->x - 1 + DIM) % DIM, s->y },
		{ s->x, (s->y - 1 + DIM) % DIM }
	};
	int sp, i;

	sp = lat[s->x][s->y];
	lat[s->x][s->y] *= -1;
	*m -= 2 * sp;
	for (i = 0; i < 4; ++i)
		if (lat[nbr[i].x][nbr[i].y] == sp && frand(0, 1) < p)
			wolff(&nbr[i], p, m);
}

int main()
{
	int i;
	long m;
	double p, t, mag, msqr, chi;
	struct Site s;
	const int steps = 1e4;

	genlat();
	m = measure();
	for (t = 2.4; t >= 2.1; t -= 0.005) {
		fprintf(stderr, "t = %lf\t", t);
		p = 1 - exp(-2.0 / t);
		mag = msqr = 0.0;
		for (i = -steps; i < steps; ++i) {
			s.x = (int) frand (0, DIM - EPS);
			s.y = (int) frand (0, DIM - EPS);
			wolff(&s, p, &m);
			if (i >= 0) {
				mag += fabs((double) m) / NSPINS;
				msqr += pow((double) m / NSPINS, 2);
			}
		}
		mag /= steps;
		msqr /= steps;
		chi = (msqr - mag * mag ) / t;
		printf("%e\t%e\t%e\n", t, mag, chi);
		fprintf(stderr, "<m> = %lf\n", mag);
	}
	return 0;
}
