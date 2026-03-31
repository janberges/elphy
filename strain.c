#include "elphy.h"

/* calculate energy of strained interatomic springs per unit cell */

double strain_energy(const struct model m) {
    struct element *k;
    double energy, ux, uy;
    div_t x, y;
    int i;

    energy = 0.0;

    for (k = m.k; k - m.k < m.nk; k++) {
        x = div(k->a, 3);
        y = div(k->b, 3);

        ux = m.tau[y.quot][x.rem] - m.tau[x.quot][x.rem];
        uy = m.tau[y.quot][y.rem] - m.tau[x.quot][y.rem];

        for (i = 0; i < 3; i++) {
            ux += m.r[k->r][i] * m.uc[i][x.rem];
            uy += m.r[k->r][i] * m.uc[i][y.rem];
        }

        energy += ux * k->c * uy;
    }

    return -0.25 * m.strain * m.strain * energy;
}

/* add change of hopping and on-site energies due to uniform strain */

void strain(double **h, const struct model m, const int nc, const int **cr) {
    struct vertex *g;
    double gu;
    div_t x;
    int i, c;

    for (g = m.g; g - m.g < m.ng; g++) {
        x = div(g->x, 3);

        gu = m.tau[x.quot][x.rem];

        for (i = 0; i < 3; i++)
            gu += m.r[g->rph][i] * m.uc[i][x.rem];

        gu *= m.strain * g->c;

        for (c = 0; c < nc; c++)
            h[m.nel * c + g->a][m.nel * cr[c][g->rel] + g->b] += gu;
    }
}
