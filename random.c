#include "elphy.h"

/* MINSTD using double instead of long long int (to avoid C99 feature) */

double minstd(void) {
    const double a = 48271.0, m = 2147483647.0;
    static double i = 1.0;

    i = fmod(i * a, m);

    return i / m;
}

void random_displacements(const int nat, double *u, double umax) {
    double rho, norm, avg;
    int i, j;

    for (i = 0; i < nat; i++) {
        rho = umax * minstd();

        norm = 0.0;

        for (j = 3 * i; j < 3 * i + 3; j++) {
            u[j] = 1.0 - 2.0 * minstd();
            norm += u[j] * u[j];
        }

        norm = sqrt(norm);

        if (norm != 0.0)
            for (j = 0; j < 3; j++)
                u[3 * i + j] *= rho / norm;
    }

    for (j = 0; j < 3; j++) {
        avg = 0.0;

        for (i = 0; i < nat; i++)
            avg += u[3 * i + j];

        avg /= nat;

        for (i = 0; i < nat; i++)
            u[3 * i + j] -= avg;
    }
}
