#include "spline.h"
#include <stddef.h>


static int find_segment(double x, int n, const double *x_nodes)
{
    if (x <= x_nodes[0])
        return 0;

    if (x >= x_nodes[n - 1])
        return n - 2;

    int left = 0;
    int right = n - 1;

    while (right - left > 1) {
        int mid = (left + right) / 2;

        if (x_nodes[mid] <= x)
            left = mid;
        else
            right = mid;
    }

    return left;
}

int spline_init(int n, const double *x, const double *f, double *coef)
{
    if (n < 4 || x == NULL || f == NULL || coef == NULL)
        return -1;

    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i - 1])
            return -2;
    }

    double *lower = new double[n];
    double *diag = new double[n];
    double *upper = new double[n];
    double *rhs = new double[n];
    double *d = new double[n];

    for (int i = 0; i < n; i++) {
        lower[i] = 0.0;
        diag[i] = 0.0;
        upper[i] = 0.0;
        rhs[i] = 0.0;
        d[i] = 0.0;
    }

    double h0 = x[1] - x[0];
    double h1 = x[2] - x[1];
    double div0 = (f[1] - f[0]) / h0;
    double div1 = (f[2] - f[1]) / h1;

    diag[0] = h1;
    upper[0] = h0 + h1;
    rhs[0] = (div0 * h1 * (2.0 * x[2] + x[1] - 3.0 * x[0]) +
              div1 * h0 * h0) /
             (h0 + h1);

    for (int i = 1; i <= n - 2; i++) {
        double hm = x[i] - x[i - 1];
        double hp = x[i + 1] - x[i];
        double divm = (f[i] - f[i - 1]) / hm;
        double divp = (f[i + 1] - f[i]) / hp;

        lower[i] = hp;
        diag[i] = 2.0 * (hm + hp);
        upper[i] = hm;
        rhs[i] = 3.0 * divm * hp + 3.0 * divp * hm;
    }

    double hl = x[n - 2] - x[n - 3];
    double hr = x[n - 1] - x[n - 2];
    double divl = (f[n - 2] - f[n - 3]) / hl;
    double divr = (f[n - 1] - f[n - 2]) / hr;

    lower[n - 1] = hl + hr;
    diag[n - 1] = hl;
    rhs[n - 1] =
        (divl * hr * hr +
         divr * hl * (3.0 * x[n - 1] - x[n - 2] - 2.0 * x[n - 3])) /
        (hl + hr);

    for (int i = 1; i < n; i++) {
        double m = lower[i] / diag[i - 1];

        diag[i] -= m * upper[i - 1];
        rhs[i] -= m * rhs[i - 1];
    }

    d[n - 1] = rhs[n - 1] / diag[n - 1];

    for (int i = n - 2; i >= 0; i--)
        d[i] = (rhs[i] - upper[i] * d[i + 1]) / diag[i];

    for (int i = 0; i < n - 1; i++) {
        double h = x[i + 1] - x[i];
        double div = (f[i + 1] - f[i]) / h;

        coef[4 * i + 0] = f[i];
        coef[4 * i + 1] = d[i];
        coef[4 * i + 2] = (3.0 * div - 2.0 * d[i] - d[i + 1]) / h;
        coef[4 * i + 3] = (d[i] + d[i + 1] - 2.0 * div) / (h * h);
    }

    delete[] lower;
    delete[] diag;
    delete[] upper;
    delete[] rhs;
    delete[] d;

    return 0;
}

double spline_compute(double x, double a, double b, int n,
                      const double *x_nodes, const double *coef)
{
    (void)a;
    (void)b;

    if (n < 2 || x_nodes == NULL || coef == NULL)
        return 0.0;

    int i = find_segment(x, n, x_nodes);
    double t = x - x_nodes[i];

    return coef[4 * i + 0] +
           coef[4 * i + 1] * t +
           coef[4 * i + 2] * t * t +
           coef[4 * i + 3] * t * t * t;
}
