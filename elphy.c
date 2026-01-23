#include <stdio.h>
#include <stdlib.h>

/* declare interface to LAPACK subroutine */

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

double **matrix(const int n);

double *eigenvalues(const int n, double **a);

struct element {
    int r, a, b;
    double c;
};

struct model {
    int nr, nb, nt;
    int (*r)[3];
    struct element *t;
};

void get_model(const char *filename, struct model *m);
void put_model(const char *filename, struct model *m);

int main(int argc, char **argv) {
    double **h, *e;
    struct model m;
    struct element *t;
    int n, nc, c1, c2, i, j, k, l;

    /* read (and write) tight-binding model */

    get_model("el.dat", &m);
    put_model("el_copy.dat", &m);

    /* get supercell size from command-line argument */

    nc = argc > 1 ? atoi(argv[1]) : 12;
    n = m.nb * nc * nc;

    /* emulate 2D variable-length array (to avoid C99 feature) */

    h = matrix(n);

    /* populate matrix using example of supercell tight-binding Hamiltonian */

    for (i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
            c1 = m.nb * (i * nc + j);

            for (t = m.t; t - m.t < m.nt; t++) {
                k = (i + m.r[t->r][0]) % nc;
                l = (j + m.r[t->r][1]) % nc;

                if (k < 0) k += nc;
                if (l < 0) l += nc;

                c2 = m.nb * (k * nc + l);

                h[c1 + t->a][c2 + t->b] = t->c;
            }
        }
    }

    /* diagonalize matrix */

    e = eigenvalues(n, h);

    /* print eigenvalues */

    for (i = 0; i < n; i++)
        printf("% .9f\n", e[i]);

    return 0;
}

double **matrix(const int n) {
    double **a, *data;
    int row;

    a = malloc(n * sizeof data);
    data = calloc(n * n, sizeof *data);

    for (row = 0; row < n; row++)
        a[row] = data + row * n;

    return a;
}

double *eigenvalues(const int n, double **a) {
    const char jobz = 'N', uplo = 'U';
    double *w, *work, lworkopt;
    int lwork, info, step;

    w = malloc(n * sizeof *w);

    for (step = 0; step < 2; step++) {
        if (step) { /* workspace allocation */
            lwork = (int) lworkopt;
            work = malloc(lwork * sizeof *work);
        } else { /* workspace query */
            lwork = -1;
            work = &lworkopt;
        }

        dsyev_(&jobz, &uplo, &n, *a, &n, w, work, &lwork, &info);
    }

    return w;
}

void get_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3];
    struct element *t;

    fp = fopen(filename, "r");

    fscanf(fp, "%d %d %d", &m->nr, &m->nb, &m->nt);

    m->r = malloc(m->nr * sizeof *r);
    m->t = malloc(m->nt * sizeof *t);

    for (r = m->r; r - m->r < m->nr; r++)
        fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2);

    for (t = m->t; t - m->t < m->nt; t++)
        fscanf(fp, "%d %d %d %lf", &t->r, &t->a, &t->b, &t->c);

    fclose(fp);
}

void put_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3];
    struct element *t;

    fp = fopen(filename, "w");

    fprintf(fp, "%d %d %d\n", m->nr, m->nb, m->nt);

    for (r = m->r; r - m->r < m->nr; r++)
        fprintf(fp, "% d % d % d\n", **r, *(*r + 1), *(*r + 2));

    for (t = m->t; t - m->t < m->nt; t++)
        fprintf(fp, "%d %d %d % g\n", t->r, t->a, t->b, t->c);

    fclose(fp);
}
