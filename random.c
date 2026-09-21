#include "elphy.h"

/* generate normally distributed random numbers (<x> = 0, <x^2> = 1) */

double box_muller() {
    const double pi = 4.0 * atan(1.0);
    static double rho, phi;
    static int havedata = 0;

    if (havedata = !havedata) {
        while (!(rho = (double) rand() / (double) RAND_MAX));
        while (!(phi = (double) rand() / (double) RAND_MAX));

        rho = sqrt(-2.0 * log(rho));
        phi *= 2.0 * pi;

        return rho * cos(phi);
    } else
        return rho * sin(phi);
}

/* generate isotropcially distributed random atomic displacements */

void random_displacements(const int nat, double *u, const double umax) {
    double norm, scale;
    int i, j;

    for (i = 0; i < nat; i++) {
        norm = 0.0;

        for (j = 3 * i; j < 3 * i + 3; j++) {
            u[j] = box_muller();
            norm += u[j] * u[j];
        }

        norm = sqrt(norm);

        if (norm != 0.0) {
            scale = umax / norm * (double) rand() / (double) RAND_MAX;

            for (j = 0; j < 3; j++)
                u[3 * i + j] *= scale;
        }
    }

    fixcom(nat, u);
}

/* fix center of mass by setting average atomic displacement to zero */

void fixcom(const int nat, double *u) {
    double avg;
    int i, j;

    for (j = 0; j < 3; j++) {
        avg = 0.0;

        for (i = 0; i < nat; i++)
            avg += u[3 * i + j];

        avg /= nat;

        for (i = 0; i < nat; i++)
            u[3 * i + j] -= avg;
    }
}
