#include "elphy.h"

/* map lattice vectors from unit cells to supercell */

int **map(const int nc, const int nr, const int (*points)[3]) {
    int **cr, c, r, i, j, k, l;

    cr = array_2d(nc * nc, nr);

    for (i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
            c = i * nc + j;

            for (r = 0; r < nr; r++) {
                k = (i + points[r][0]) % nc;
                l = (j + points[r][1]) % nc;

                if (k < 0) k += nc;
                if (l < 0) l += nc;

                cr[c][r] = k * nc + l;
            }
        }
    }

    return cr;
}

/* populate matrix using example of supercell tight-binding Hamiltonian */

void supercell(double **h, const struct model m, const int nc, const int **cr) {
    struct element *t;
    int c;

    for (c = 0; c < nc * nc; c++)
        for (t = m.t; t - m.t < m.nt; t++)
            h[m.nb * c + t->a][m.nb * cr[c][t->r] + t->b] = t->c;
}

/* add linear electron-lattice coupling to supercell Hamiltonian */

void perturbation(double **h, const struct coupling m, const double *u,
    const int nc, const int **cr) {

    struct vertex *g;
    int c;

    for (c = 0; c < nc * nc; c++)
        for (g = m.g; g - m.g < m.ng; g++)
            h[m.nel * c + g->a][m.nel * cr[c][g->rel] + g->b]
                += u[m.nph * cr[c][g->rph] + g->x] * g->c;
}

/* calculate Jacobian via Hellmann-Feynman theorem */

double *jacobian(const double **h, const struct coupling m, const double *occ,
    const int nc, const int **cr) {

    struct vertex *g;
    int c, n, i0, iel, iph;
    double *f;

    f = calloc(m.nph * nc * nc, sizeof *f);

    for (c = 0; c < nc * nc; c++)
        for (g = m.g; g - m.g < m.ng; g++) {
            i0 = m.nel * c + g->a;
            iel = m.nel * cr[c][g->rel] + g->b;
            iph = m.nph * cr[c][g->rph] + g->x;

            for (n = 0; n < m.nel * nc * nc; n++)
                f[iph] += g->c * h[n][i0] * occ[n] * h[n][iel];
        }

    return f;
}
