#ifndef HERMITE_H
#define HERMITE_H

int hermite_init(int n, const double *x, const double *f,
                 const double *df, double *coef);

double hermite_compute(double x, double a, double b, int n,
                       const double *x_nodes, const double *coef);

#endif
