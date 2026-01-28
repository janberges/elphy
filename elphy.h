#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

double **matrix(const int n);

int **array_2d(const int rows, const int cols);

int eigenvalues(const int n, double **a, double *w);

struct element {
    int r, a, b;
    double c;
};

struct vertex {
    int rph, x, rel, a, b;
    double c;
};

struct model {
    double kt, n;
    int sc[3][3];
    double uc[3][3];
    double (*tau)[3];
    char (*typ)[3];
    int nel, nph, nat, nr, nt, nk, ng;
    int (*r)[3];
    struct element *t;
    struct element *k;
    struct vertex *g;
};

void get_model(const char *filename, struct model *m);

void get_displ(const char *filename, const int nat,
    double uc[3][3], char (*typ)[3], double (*tau)[3], double *u);

void put_displ(const char *filename, const int nat,
    double uc[3][3], char (*typ)[3], double (*tau)[3], double *u);

int dot(const int a[3], const int b[3]);

void cross(const int a[3], const int b[3], int *c);

int map(const struct model m, int ***cr, int ***cells);

void repeat(const struct model m, const int nc, int **cells,
    double uc[3][3], char (*typ)[3], double (*tau)[3]);

void supercell(double **a, const int nb, const int nl, const struct element *l,
    const int nc, int **cr);

void perturbation(double **h0, double **h, const struct model m,
    const double *u, const int nc, int **cr);

void jacobian(double **h, const struct model m, const double *occ,
    double *f, const int nc, int **cr);

double fermi(const double x);
double dirac(const double x);

double fermi_level(const int ne, const double n,
    const double *e, const double kt, double mu);

double free_energy(const int ne, const double n,
    const double *e, const double kt, const double mu);

void occupations(const int ne, double *f,
    const double *e, const double kt, const double mu);
