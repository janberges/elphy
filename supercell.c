#include "elphy.h"

static int dot(const int a[3], const int b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross(const int a[3], const int b[3], int *c) {
    int i;
    for (i = 0; i < 3; i++)
        c[i] = a[(i + 1) % 3] * b[(i + 2) % 3]
             - a[(i + 2) % 3] * b[(i + 1) % 3];
}

/* find unit-cell vectors within supercell and add lattice vectors */

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

void repeat(double uc[3][3], char (*typ)[3], double (*tau)[3], double (*fdc)[3],
    const struct model m, const int nc, int **cells) {

    int i, j, k, c, ci;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            uc[i][j] = 0.0;
            for (k = 0; k < 3; k++)
                uc[i][j] += m.sc[i][k] * m.uc[k][j];
        }

    for (c = 0; c < nc; c++)
        for (i = 0; i < m.nat; i++) {
            ci = m.nat * c + i;
            for (j = 0; j < 3; j++) {
                typ[ci][j] = m.typ[i][j];
                tau[ci][j] = m.tau[i][j];
                fdc[ci][j] = m.fdc[i][j];
                for (k = 0; k < 3; k++)
                    tau[ci][j] += cells[c][k] * m.uc[k][j];
            }
        }
}

/* populate supercell tight-binding Hamiltonian or force-constants matrix */

void populate(double **a, const int nb, const int nl, const struct element *l,
    const int nc, int **cr) {

    const struct element *m;
    int c;

    for (c = 0; c < nc; c++)
        for (m = l; m - l < nl; m++)
            a[nb * c + m->a][nb * cr[c][m->r] + m->b] += m->c;
}

/* add linear electron-lattice coupling to supercell Hamiltonian */

void perturb(double **h, const double *u, const struct model m,
    const int nc, int **cr) {

    struct vertex *g;
    int c;

    for (c = 0; c < nc; c++)
        for (g = m.g; g - m.g < m.ng; g++)
            h[m.nel * c + g->a][m.nel * cr[c][g->rel] + g->b]
                += u[m.nph * cr[c][g->rph] + g->x] * g->c;
}

/* add electronic contribution to supercell forces */

void add_forces(double *forces, double **occ, const struct model m,
    const int nc, int **cr) {

    struct vertex *g;
    int c;

    for (c = 0; c < nc; c++)
        for (g = m.g; g - m.g < m.ng; g++)
            forces[m.nph * cr[c][g->rph] + g->x]
                -= occ[m.nel * c + g->a][m.nel * cr[c][g->rel] + g->b] * g->c;
}
