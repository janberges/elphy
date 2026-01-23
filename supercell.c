#include "elphy.h"

/* populate matrix using example of supercell tight-binding Hamiltonian */

void supercell(double **h, struct model m, int nc) {
    struct element *t;
    int c0, cr, i, j, k, l;

    for (i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
            c0 = m.nb * (i * nc + j);

            for (t = m.t; t - m.t < m.nt; t++) {
                k = (i + m.r[t->r][0]) % nc;
                l = (j + m.r[t->r][1]) % nc;

                if (k < 0) k += nc;
                if (l < 0) l += nc;

                cr = m.nb * (k * nc + l);

                h[c0 + t->a][cr + t->b] = t->c;
            }
        }
    }
}

/* add linear electron-lattice coupling to supercell Hamiltonian */

void perturbation(double **h, struct coupling m, double *u, int nc) {
    struct vertex *g;
    int c0, crel, crph, i, j, k, l, p, q;

    for (i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
            c0 = m.nel * (i * nc + j);

            for (g = m.g; g - m.g < m.ng; g++) {
                k = (i + m.r[g->rel][0]) % nc;
                l = (j + m.r[g->rel][1]) % nc;

                p = (i + m.r[g->rph][0]) % nc;
                q = (j + m.r[g->rph][1]) % nc;

                if (k < 0) k += nc;
                if (l < 0) l += nc;

                if (p < 0) p += nc;
                if (q < 0) q += nc;

                crel = m.nel * (k * nc + l);
                crph = m.nph * (p * nc + q);

                h[c0 + g->a][crel + g->b] += u[crph + g->x] * g->c;
            }
        }
    }
}

/* calculate Jacobian via Hellmann–Feynman theorem */

double *jacobian(double **h, struct coupling m, double *occ, int nc) {
    struct vertex *g;
    int c0, crel, crph, i, j, k, l, p, q, n;
    double *f;

    f = calloc(m.nph * nc * nc, sizeof *f);

    for (i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
            c0 = m.nel * (i * nc + j);

            for (g = m.g; g - m.g < m.ng; g++) {
                k = (i + m.r[g->rel][0]) % nc;
                l = (j + m.r[g->rel][1]) % nc;

                p = (i + m.r[g->rph][0]) % nc;
                q = (j + m.r[g->rph][1]) % nc;

                if (k < 0) k += nc;
                if (l < 0) l += nc;

                if (p < 0) p += nc;
                if (q < 0) q += nc;

                crel = m.nel * (k * nc + l);
                crph = m.nph * (p * nc + q);

                for (n = 0; n < m.nel * nc * nc; n++)
                    f[crph + g->x] += g->c
                        * h[n][c0 + g->a] * occ[n] * h[n][crel + g->b];
            }
        }
    }

    return f;
}
