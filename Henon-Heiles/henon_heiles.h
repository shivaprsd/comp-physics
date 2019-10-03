#ifndef HENON_HEILES_H
#define HENON_HEILES_H

double hh(double x, double y);
void hh_evol(double t, double v[4], double vdot[4]);
void hh_tan(double t, double w[20], double wdot[20]);
double hh_equi(double y, double u);
double px(double e, double x, double y, double py);
double en(double *x, int n);

void hh_map(double xin[2], double xout[2], double a);
void hh_map_tan(double win[6], double wout[6], double a);

void plot_equi();
void plot_map(double *x, double a, int n);
void phase_space(double *x, double energy, int n);

#endif
