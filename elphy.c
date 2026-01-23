#include <math.h>
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

void supercell(double **h, struct model m, int nc);

double fermi(double x);
double dirac(double x);

double fermi_level(const int ne, const double n,
    const double *e, const double kt, double mu);

double free_energy(const int ne, const double n,
    const double *e, const double kt, const double mu);

int main(int argc, char **argv) {
    const double kt = 0.025;
    double **h, *e, nel = 0.25, mu = 0.0;
    struct model m;
    int n, nc;

    /* read (and write) tight-binding model */

    get_model("el.dat", &m);
    put_model("el_copy.dat", &m);

    /* get supercell size from command-line argument */

    nc = argc > 1 ? atoi(argv[1]) : 12;
    n = m.nb * nc * nc;
    nel *= n;

    /* emulate 2D variable-length array (to avoid C99 feature) */

    h = matrix(n);

    /* populate matrix using example of supercell tight-binding Hamiltonian */

    supercell(h, m, nc);

    /* diagonalize matrix */

    e = eigenvalues(n, h);

    /* determine chemical potential */

    mu = fermi_level(n, nel, e, kt, mu);

    /* determine free energy */

    printf("%.9f\n", 2.0 * free_energy(n, nel, e, kt, mu));

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

void supercell(double **h, struct model m, int nc) {
    struct element *t;
    int c1, c2, i, j, k, l;

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
}

double fermi(double x) {
    return 0.5 - 0.5 * tanh(0.5 * x);
}

double dirac(double x) {
    return 1.0 / (2.0 * cosh(x) + 2.0);
}

double fermi_level(const int ne, const double n,
    const double *e, const double kt, double mu) {

    const double eps = 1e-10, tol = 1e-5;
    const double f0 = fermi(0.0);
    const double d0 = dirac(0.0) / kt;
    double x, f, w, sum_f, sum_w, sum_e_w;
    int i;

    for (;;) {
        sum_f = sum_w = sum_e_w = 0.0;

        for (i = 0; i < ne; i++) {
            x = e[i] - mu;
            f = fermi(x / kt);

            w = fabs(x) > eps ? (f0 - f) / x : d0;

            sum_f += f;
            sum_w += w;
            sum_e_w += e[i] * w;
        }

        if (fabs(sum_f- n) < tol)
            return mu;

        mu = (n - ne * f0 + sum_e_w) / sum_w;
    }
}

double free_energy(const int ne, const double n,
    const double *e, const double kt, const double mu) {

    double grand = 0.0;
    int i;

    for (i = 0; i < ne; i++)
        grand -= log(exp((mu - e[i]) / kt) + 1.0);

    grand *= kt;

    return grand + n * mu;
}
