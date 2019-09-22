#include <stdio.h>
#include "henon_heiles.c"

/* Program to calculate the relative area covered by invariant curves in the
 * phase space section (x = 0, px > 0) for different energies, by Monte-Carlo
 * method (using 100 random points). Also calculates the Lapunov Exponents for
 * the first three trajectories for each energy */
int main()
{
	int pc, invc, i, j;
	long t;
	double lda, mu, muc = 1e-4, h = 1e-3;
	double x1[4], x1new[4], x2[4];
	double einv[] = { 50.0, 25.0, 16.0, 12.0, 10.5, 10.0, 9.5, 9.0, 8.5,
			8.0, 7.0, 6.5, 6.0 };
	const int np = 100;

	srand(time(NULL));	/* initialise rand() with system time */
	printf("#Energy\t\tlambda_1\tlambda_2\tlambda_3\tRelative area\n");
	for (j = 0; j < LEN(einv); ++j) {
		printf("%lf\t", 1 / einv[j]);
		for (invc = i = 0; i < np; ++i) {
			x1[0] = 0.0;
			/* Generate random points until we get a valid one */
			do {
				x1[1] = frand(-0.5, 1.0);
				x1[3] = frand(-0.5, 0.5);
				x1[2] = px(1.0 / einv[j], x1[0], x1[1], x1[3]);
			} while (x1[2] < 0);

			copy_vec(x1, x2, 4);
			x2[1] += 1e-7;	/* Make x2 slightly separated from x1 */
			mu = 0.0;
			t = pc = 0;
			/* Calculate mu for 25 transformations of x */
			while (pc < 25) {
				rk4(h * t, x1, 4, hh_evol, h, x1new);
				rk4(h * t++, x2, 4, hh_evol, h, x2);
				/* Again, check if we crossed the x = 0 plane */
				if (x1new[0] * x1[0] < 0 && x1[2] > 0.0) {
					mu += pow(dist(x1, x2, 4), 2);
					++pc;
				}
				copy_vec(x1new, x1, 4);
			}
			if (mu < muc)
				++invc;
			/* Print lambda for only three trajectories */
			if (i < 3) {
				lda = log(dist(x1, x2, 4) * 1e7) / (h * t);
				printf("%e\t", lda);
			}
		}
		/* Print ratio of invariant trajectories to the total */
		printf("%lf\n", invc / (float) np);
	}
	return 0;
}

