#include <stdio.h>
#include "henon_heiles.c"

/* Program to closely reproduce the map obtained by Hénon and Heiles (fig. 8) */
int main()
{
	int i;
	/* Initial points */
	double fig8[][2] = {
		{ 0.0, 0.2 }, { 0.0, 0.4 }, { 0.0, 0.48 }, { 0.0, 0.5 },
		{ 0.0, 0.7 }, { 0.8, 0.0 }, { -0.3, -0.3 }, { 0.3, 0.3 },
		{ 0.0, 0.7 }, { 0.8, 0.0 }, { -0.3, -0.3 }, { 0.3, 0.3 },
		{ 0.0, 0.615 }, { 0.615, 0.0 }, { 0.0, 0.0 }
	};

	for (i = 0; i < LEN(fig8); ++i) {
		plot_map(fig8[i], 1.6, 2000);
		printf("\n\n");
	}
	return 0;
}
