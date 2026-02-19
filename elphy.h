#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

void dsymv_(const char *uplo, const int *n, const double *alpha,
    const double *a, const int *lda, const double *x, const int *incx,
    const double *beta, double *y, const int *incy);

double ddot_(const int *n, const double *dx, const int *incx,
    const double *dy, const int *incy);

void daxpy_(const int *n, const double *da, const double *dx, const int *incx,
    double *dy, const int *incy);

void dscal_(const int *n, const double *da, double *dx, const int *incx);

void dsyrk_(const char *uplo, const char *trans, const int *n,
    const int *k, const double *alpha, const double *a, const int *lda,
    const double *beta, double *c, const int *ldc);

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
    int nel, nspin;
    int sc[3][3];
    double uc[3][3];
    int nph, nat;
    char (*typ)[3];
    double (*tau)[3];
    double (*fdc)[3];
    int nr, (*r)[3];
    int nt, nk, ng;
    struct element *t;
    struct element *k;
    struct vertex *g;
};

double step(double **h, double **h0, double *e, double **occ,
    double **c, const double *u, double *forces, const double *forces0,
    const struct model m, const int nc, int **cr,
    const int lwork, double *work);

void driver(char *host, double **h, double **h0, double *e, double **occ,
    double **c, double *u, double *forces, const double *forces0,
    double (*tau)[3], const struct model m, const int nc, int **cr,
    const int lwork, double *work);

void error(char *msg, ...);

void get_model(const char *filename, struct model *m);

int get_displ(const int nat, char (*typ)[3], double (*tau)[3], double *u);

void put_displ(const char *filename, const int nat, double uc[3][3],
    char (*typ)[3], double (*tau)[3], double *u);

void put_force(const int nat, char (*typ)[3],
    const double energy, const double *forces);

double **matrix(const int n);

int **array_2d(const int rows, const int cols);

void random_displacements(const int nat, double *u, double umax);

int open_inet_socket(const char *host, const char *port);

int open_unix_socket(const char *host, const char *prefix);

void sread(const int sfd, void *data, const int len);

void swrite(const int sfd, const void *data, const int len);

int map(const struct model m, int ***cr, int ***cells);

void repeat(double uc[3][3], char (*typ)[3], double (*tau)[3], double (*fdc)[3],
    const struct model m, const int nc, int **cells);

void populate(double **a, const int nb, const int nl, const struct element *l,
    const int nc, int **cr);

void perturb(double **h, const double *u, const struct model m,
    const int nc, int **cr);

void add_forces(double *forces, double **occ, const struct model m,
    const int nc, int **cr);

double fermi_level(double n, const int ne, const double *e, const double kt,
    double mu);

double grand_potential(const int ne, const double *e, const double kt,
    const double mu);

void occupations(const int ne, const double *e, const double kt,
    const double mu, const int nspin, double **psi, double **occ);
