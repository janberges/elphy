#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

struct vertex {
    int rph, x, rel, a, b;
    double c;
};

struct coupling {
    int nr, nel, nph, ng;
    int (*r)[3];
    struct vertex *g;
};

void get_model(const char *filename, struct model *m);
void put_model(const char *filename, struct model *m);
void get_coupl(const char *filename, struct coupling *m);
void get_displ(const char *filename, const int nx, double *u);

void supercell(double **h, struct model m, int nc);

void perturbation(double **h, struct coupling m, double *u, int nc);

double fermi(double x);
double dirac(double x);

double fermi_level(const int ne, const double n,
    const double *e, const double kt, double mu);

double free_energy(const int ne, const double n,
    const double *e, const double kt, const double mu);

int main(int argc, char **argv) {
    const double kt = 0.0019;
    double **h, **c, *e, *u, energy, nel = 0.25, mu = 0.0;
    struct model el, ph;
    struct coupling elph;
    int n, nx, nc, i, j;

    get_model("el.dat", &el);
    get_model("ph.dat", &ph);
    get_coupl("elph.dat", &elph);

    /* put_model("el_copy.dat", &el); */

    nc = argc > 1 ? atoi(argv[1]) : 12;
    nx = ph.nb * nc * nc;
    n = el.nb * nc * nc;
    nel *= n;

    h = matrix(n);
    c = matrix(nx);

    supercell(h, el, nc);
    supercell(c, ph, nc);

    u = malloc(nx * sizeof *u);
    get_displ("u.dat", nx, u);

    perturbation(h, elph, u, nc);

    e = eigenvalues(n, h);

    mu = fermi_level(n, nel, e, kt, mu);

    energy = 2.0 * free_energy(n, nel, e, kt, mu);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    printf("%.9f\n", energy);

    return 0;
}

/* emulate 2D variable-length array (to avoid C99 feature) */

double **matrix(const int n) {
    double **a, *data;
    int row;

    a = malloc(n * sizeof data);
    data = calloc(n * n, sizeof *data);

    for (row = 0; row < n; row++)
        a[row] = data + row * n;

    return a;
}

/* diagonalize matrix using LAPACK subroutine */

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

/* read tight-binding or mass-spring model */

void get_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3];
    struct element *t;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run data.py first.\n", filename);
        exit(1);
    }

    fscanf(fp, "%d %d %d", &m->nb, &m->nr, &m->nt);

    m->r = malloc(m->nr * sizeof *r);
    m->t = malloc(m->nt * sizeof *t);

    for (r = m->r; r - m->r < m->nr; r++)
        fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2);

    for (t = m->t; t - m->t < m->nt; t++)
        fscanf(fp, "%d %d %d %lf", &t->r, &t->a, &t->b, &t->c);

    fclose(fp);
}

/* write tight-binding or mass-spring model */

void put_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3];
    struct element *t;

    fp = fopen(filename, "w");

    fprintf(fp, "%d %d %d\n", m->nb, m->nr, m->nt);

    for (r = m->r; r - m->r < m->nr; r++)
        fprintf(fp, "% d % d % d\n", **r, *(*r + 1), *(*r + 2));

    for (t = m->t; t - m->t < m->nt; t++)
        fprintf(fp, "%d %d %d % .9f\n", t->r, t->a, t->b, t->c);

    fclose(fp);
}

/* read coupling model */

void get_coupl(const char *filename, struct coupling *m) {
    FILE *fp;
    int (*r)[3];
    struct vertex *g;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run data.py first.\n", filename);
        exit(1);
    }

    fscanf(fp, "%d %d %d %d", &m->nph, &m->nel, &m->nr, &m->ng);

    m->r = malloc(m->nr * sizeof *r);
    m->g = malloc(m->ng * sizeof *g);

    for (r = m->r; r - m->r < m->nr; r++)
        fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2);

    for (g = m->g; g - m->g < m->ng; g++)
        fscanf(fp, "%d %d %d %d %d %lf",
            &g->rph, &g->x, &g->rel, &g->a, &g->b, &g->c);

    fclose(fp);
}

/* read atomic displacements */

void get_displ(const char *filename, const int nx, double *u) {
    FILE *fp;
    int i;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run test.py first.\n", filename);
        exit(1);
    }

    for (i = 0; i < nx; i++)
        fscanf(fp, "%lf", u++);

    fclose(fp);
}

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
