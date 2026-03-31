#include "elphy.h"

/* calculate lattice energy due to uniform strain */

double strain_energy(const struct model m) {
    struct element *k;
    double energy, ux, uy;
    div_t x, y;
    int i;

    energy = 0.0;

    for (k = m.k; k - m.k < m.nk; k++) {
        x = div(k->a, 3);
        y = div(k->b, 3);

        ux = m.tau[x.quot][x.rem] - m.tau[y.quot][x.rem];
        uy = m.tau[y.quot][y.rem] - m.tau[x.quot][y.rem];

        for (i = 0; i < 3; i++) {
            ux -= m.r[k->r][i] * m.uc[i][x.rem];
            uy += m.r[k->r][i] * m.uc[i][y.rem];
        }

        energy += ux * k->c * uy;
    }

    return 0.25 * m.strain * m.strain * energy;
}
