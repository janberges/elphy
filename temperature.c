#include "elphy.h"

static double fermi(const double x) {
    return 0.5 - 0.5 * tanh(0.5 * x);
}

static double dirac(const double x) {
    return 1.0 / (2.0 * cosh(x) + 2.0);
}

double fermi_level(double n, const int ne, const double *e, const double kt,
    double mu) {

    const double eps = 1e-10, tol = 1e-5;
    const double f0 = fermi(0.0);
    const double d0 = dirac(0.0) / kt;
    double x, f, w, sum_f, sum_w, sum_e_w;
    int i;

    for (;;) {
        sum_f = sum_w = sum_e_w = 0.0;

        for (i = 0; i < ne; i++) {
            x = e[i] - mu;
            f = fermi(x / kt);

            w = fabs(x) > eps ? (f0 - f) / x : d0;

            sum_f += f;
            sum_w += w;
            sum_e_w += e[i] * w;
        }

        if (fabs(sum_f - n) < tol)
            return mu;

        mu = (n - ne * f0 + sum_e_w) / sum_w;
    }
}

double grand_potential(const int ne, const double *e, const double kt,
    const double mu) {

    double grand = 0.0, x;
    int i;

    for (i = 0; i < ne; i++) {
        x = (mu - e[i]) / kt;
        grand -= x > 709.0 ? x : log(exp(x) + 1.0);
    }

    return kt * grand;
}

void occupations(const int ne, const double *e, const double kt,
    const double mu, const int nspin, double **psi, double **occ) {

    const int inc = 1;
    const double zero = 0.0, one = 1.0;
    double tmp;
    int n, a, b;

    for (n = 0; n < ne; n++) {
        tmp = sqrt(nspin * fermi((e[n] - mu) / kt));
        dscal_(&ne, &tmp, psi[n], &inc);
    }

    dsyrk_("U", "N", &ne, &ne, &one, *psi, &ne, &zero, *occ, &ne);

    for (a = 0; a < ne; a++)
        for (b = 0; b < a; b++)
            occ[b][a] = occ[a][b];
}
