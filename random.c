#include "elphy.h"

void random_displacements(const int nat, double *u, double umax) {
    double rho, norm, avg;
    int i, j;

    for (i = 0; i < nat; i++) {
        rho = umax * (double) rand() / (double) RAND_MAX;

        norm = 0.0;

        for (j = 3 * i; j < 3 * i + 3; j++) {
            u[j] = 1.0 - 2.0 * (double) rand() / (double) RAND_MAX;
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
