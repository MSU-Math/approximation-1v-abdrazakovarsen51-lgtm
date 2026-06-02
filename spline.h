#ifndef SPLINE_H
#define SPLINE_H

int spline_init(int n, const double *x, const double *f, double *coef);

double spline_compute(double x, double a, double b, int n,
                      const double *x_nodes, const double *coef);

#endif
