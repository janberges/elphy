#include "elphy.h"

int dot(const int a[3], const int b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void cross(const int a[3], const int b[3], int *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/* map lattice vectors from unit cells to supercell */

int map(const struct model m, int ***cr, int ***cells) {
    int b[3][3], cell[3], n[3], nc, i, j, c, r, tmp;
    int lower[3] = {0, 0, 0};
    int upper[3] = {0, 0, 0};

    cross(m.sc[1], m.sc[2], b[0]);
    cross(m.sc[2], m.sc[0], b[1]);
    cross(m.sc[0], m.sc[1], b[2]);

    nc = dot(m.sc[0], b[0]);

    if (nc < 0) {
        nc *= -1;
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                b[i][j] *= -1;
    }

    for (n[0] = 0; n[0] < 2; n[0]++)
        for (n[1] = 0; n[1] < 2; n[1]++)
            for (n[2] = 0; n[2] < 2; n[2]++)
                for (j = 0; j < 3; j++) {
                    tmp = 0;
                    for (i = 0; i < 3; i++)
                        tmp += n[i] * m.sc[i][j];
                    if (tmp < lower[j])
                        lower[j] = tmp;
                    if (tmp > upper[j])
                        upper[j] = tmp;
                }

    *cells = array_2d(nc, 3);

    c = 0;
    for (n[0] = lower[0]; n[0] < upper[0]; n[0]++)
        for (n[1] = lower[1]; n[1] < upper[1]; n[1]++)
            for (n[2] = lower[2]; n[2] < upper[2]; n[2]++) {
                for (i = 0; i < 3; i++)
                    if ((tmp = dot(n, b[i])) < 0 || tmp >= nc)
                        break;

                if (i == 3) {
                    for (j = 0; j < 3; j++)
                        (*cells)[c][j] = n[j];
                    c++;
                }
            }

    if (c != nc)
        error("Problem with supercell.");

    *cr = array_2d(nc, m.nr);

    for (c = 0; c < nc; c++) {
        for (r = 0; r < m.nr; r++) {
            for (j = 0; j < 3; j++)
                cell[j] = (*cells)[c][j] + m.r[r][j];

            for (i = 0; i < 3; i++) {
                n[i] = dot(cell, b[i]) % nc;
                if (n[i] < 0)
                    n[i] += nc;
            }

            for (j = 0; j < 3; j++) {
                cell[j] = 0;
                for (i = 0; i < 3; i++)
                    cell[j] += n[i] * m.sc[i][j];
                cell[j] /= nc;
            }

            for (tmp = 0; tmp < nc; tmp++) {
                for (j = 0; j < 3; j++)
                    if (cell[j] != (*cells)[tmp][j])
                        break;

                if (j == 3) {
                    (*cr)[c][r] = tmp;
                    break;
                }
            }
        }
    }

    return nc;
}

/* determine dimensions and basis atoms of supercell */

void repeat(const struct model m, const int nc, int **cells,
    double uc[3][3], char (*typ)[3], double (*tau)[3], double (*fdc)[3]) {

    int i, j, k, c;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            uc[i][j] = 0.0;
            for (k = 0; k < 3; k++)
                uc[i][j] += m.sc[i][k] * m.uc[k][j];
        }

    for (c = 0; c < nc; c++)
        for (i = 0; i < m.nat; i++)
            for (j = 0; j < 3; j++) {
                typ[m.nat * c + i][j] = m.typ[i][j];
                tau[m.nat * c + i][j] = m.tau[i][j];
                fdc[m.nat * c + i][j] = m.fdc[i][j];
                for (k = 0; k < 3; k++)
                    tau[m.nat * c + i][j] += cells[c][k] * m.uc[k][j];
            }
}

/* populate matrix using example of supercell tight-binding Hamiltonian */

void supercell(double **a, const int nb, const int nl, const struct element *l,
    const int nc, int **cr) {

    const struct element *m;
    int c;

    for (c = 0; c < nc; c++)
        for (m = l; m - l < nl; m++)
            a[nb * c + m->a][nb * cr[c][m->r] + m->b] = m->c;
}

/* add linear electron-lattice coupling to supercell Hamiltonian */

void perturbation(double **h0, double **h, const struct model m,
    const double *u, const int nc, int **cr) {

    struct vertex *g;
    int c;

    memcpy(*h, *h0, m.nel * m.nel * nc * nc * sizeof **h);

    for (c = 0; c < nc; c++)
        for (g = m.g; g - m.g < m.ng; g++)
            h[m.nel * c + g->a][m.nel * cr[c][g->rel] + g->b]
                += u[m.nph * cr[c][g->rph] + g->x] * g->c;
}

/* calculate Jacobian via Hellmann-Feynman theorem */

void compute_forces(double **h, const struct model m, const double *occ,
    double *forces, const double *forces0, const int nc, int **cr) {

    struct vertex *g;
    int c, n, i0, iel, iph;

    memcpy(forces, forces0, nc * m.nph * sizeof *forces);

    for (c = 0; c < nc; c++)
        for (g = m.g; g - m.g < m.ng; g++) {
            i0 = m.nel * c + g->a;
            iel = m.nel * cr[c][g->rel] + g->b;
            iph = m.nph * cr[c][g->rph] + g->x;

            for (n = 0; n < m.nel * nc; n++)
                forces[iph] -= g->c * h[n][i0] * occ[n] * h[n][iel];
        }
}
