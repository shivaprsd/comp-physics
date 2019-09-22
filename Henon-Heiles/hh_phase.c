#include <stdio.h>
#include "henon_heiles.c"

/* Program to closely reproduce the phase space diagrams obtained by Hénon and
 * Heiles (fig. 4 to 6) for three different energies */
int main()
{
	int i;
	double egy, x[4];
	/* Initial points for the three figures */
	double fig4[][2] = {
		{ -0.365, 0.0 }, { -0.35, 0.0 }, { -0.3614, 0.0 },
		{ -0.29, 0.0 }, { -0.125, 0.0 }, { -0.12, 0.0 }, { 0.0, 0.0 },
		{ 0.2, 0.0 }, { 0.0, 0.24 }, { 0.0, -0.24 }
	};
	double fig5[][2] = {
		{ -0.41, 0.0 }, { -0.3614, 0.0 }, { -0.1, 0.0 }, { 0.0, 0.0 },
		{ 0.1, 0.0 }, { 0.2, 0.0 }, { 0.562, 0.0 }, { 0.0, 0.5 - EPS },
		{ 0.0, -0.225 }, { 0.0, -0.2 }, { 0.0, 0.2 }, { 0.0, 0.225 }
	};
	double fig6[][2] = {
		{ 0.0, 0.0 }, { 0.5, 0.0 }, { 1.0 - EPS, 0.0 }, { 0.3, 0.0 },
		{ -0.2, -0.44 }, { -0.2, 0.44 }
	};

	/* Demonstrating for energy = 1 / 12.0 (fig. 4) */
	/* Change values and variables correspondingly for the remaining two. */
	egy = 1 / 12.0;
	for (i = 0; i < LEN(fig4); ++i) {
		x[0] = 0.0;
		x[1] = fig4[i][0];
		x[3] = fig4[i][1];
		x[2] = px(egy, x[0], x[1], x[3]);
		phase_space(x, egy, 1000);
		printf("\n\n");
	}
	return 0;
}

