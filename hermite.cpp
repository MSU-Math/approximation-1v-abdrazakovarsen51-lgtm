#include "hermite.h"

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

int hermite_init(int n, const double *x, const double *f,
                 const double *df, double *coef)
{
    if (n < 4 || x == nullptr || f == nullptr || df == nullptr ||
        coef == nullptr)
        return -1;

    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i - 1])
            return -2;
    }

    double *d = new double[n];

    for (int i = 0; i < n; i++)
        d[i] = df[i];

    double h0 = x[1] - x[0];
    double h1 = x[2] - x[1];
    double div0 = (f[1] - f[0]) / h0;
    double div1 = (f[2] - f[1]) / h1;

    d[0] = h0 * h0 / (h1 * h1) * (d[1] + d[2] - 2.0 * div1) -
           d[1] + 2.0 * div0;

    double hl = x[n - 2] - x[n - 3];
    double hr = x[n - 1] - x[n - 2];
    double divl = (f[n - 2] - f[n - 3]) / hl;
    double divr = (f[n - 1] - f[n - 2]) / hr;

    d[n - 1] = hr * hr / (hl * hl) *
                   (d[n - 3] + d[n - 2] - 2.0 * divl) -
               d[n - 2] + 2.0 * divr;

    for (int i = 0; i < n - 1; i++) {
        double h = x[i + 1] - x[i];
        double div = (f[i + 1] - f[i]) / h;

        coef[4 * i + 0] = f[i];
        coef[4 * i + 1] = d[i];
        coef[4 * i + 2] = (3.0 * div - 2.0 * d[i] - d[i + 1]) / h;
        coef[4 * i + 3] = (d[i] + d[i + 1] - 2.0 * div) / (h * h);
    }

    delete[] d;

    return 0;
}

double hermite_compute(double x, double a, double b, int n,
                       const double *x_nodes, const double *coef)
{
    if (n < 2 || x_nodes == nullptr || coef == nullptr)
        return 0.0;

    if (x < a)
        x = a;

    if (x > b)
        x = b;

    int i = find_segment(x, n, x_nodes);
    double t = x - x_nodes[i];

    return coef[4 * i + 0] +
           coef[4 * i + 1] * t +
           coef[4 * i + 2] * t * t +
           coef[4 * i + 3] * t * t * t;
}
