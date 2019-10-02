/* Utilities for general purpose programming */

#ifndef UTIL_H
#define UTIL_H

#define LEN(arr) (sizeof(arr) / sizeof(arr[0]))
#define EPS 1e-5
#define DECK_SIZE 32    /* size of the random number shuffle table */

/* Function pointer to f(x, y1, y2...) */
typedef void (*deriv_func)(double, double [], double []);
enum bool { false, true };

double *vector(int n);
double dot(double *u, double *v, int n);
void scale(double *v, int n, double k);
double *proj(double *v, double *u, int n);
void sub(double *v, double *u, int n);
double norm(double *v, int n);
void print_vec(double *v, int n);
void copy_vec(double *s, double *t, int n);

double dist(double *v1, double *v2, int n);
double frand(double lb, double ub);
enum bool eq(double x, double y);

void rk4(double x, double y[], int n, deriv_func f, double h, double yout[]);
void ortho(double **u, int m, int n, double *a);

#endif
