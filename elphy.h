#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

extern void dsymv_(const char *uplo, const int *n, const double *alpha,
    const double *a, const int *lda, const double *x, const int *incx,
    const double *beta, double *y, const int *incy);

extern double ddot_(const int *n, const double *dx, const int *incx,
    const double *dy, const int *incy);

extern void daxpy_(const int *n, const double *da, const double *dx,
    const int *incx, double *dy, const int *incy);

extern void dscal_(const int *n, const double *da, double *dx, const int *incx);

extern void dsyrk_(const char *uplo, const char *trans, const int *n,
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
    int sc[3][3];
    double uc[3][3];
    double (*tau)[3];
    double (*fdc)[3];
    char (*typ)[3];
    int nspin, nel, nph, nat, nr, nt, nk, ng;
    int (*r)[3];
    struct element *t;
    struct element *k;
    struct vertex *g;
};

double step(double **h0, double **h, double **c, const struct model m,
    const double *u, double *e, double **occ, double *forces,
    const double *forces0, const int nc, int **cr, const int lwork,
    double *work);

void driver(double **h0, double **h, double **c, const struct model m,
    double *u, double *e, double **occ, double *forces, const double *forces0,
    double (*tau)[3], const int nc, int **cr, const int lwork, double *work,
    char *host);

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

int open_inet_socket(const char *host, int port);

int open_unix_socket(const char *host, const char *prefix);

void sread(const int sfd, void *data, const int len);

void swrite(const int sfd, const void *data, const int len);

int dot(const int a[3], const int b[3]);

void cross(const int a[3], const int b[3], int *c);

int map(const struct model m, int ***cr, int ***cells);

void repeat(const struct model m, const int nc, int **cells,
    double uc[3][3], char (*typ)[3], double (*tau)[3], double (*fdc)[3]);

void supercell(double **a, const int nb, const int nl, const struct element *l,
    const int nc, int **cr);

void perturbation(double **h0, double **h, const struct model m,
    const double *u, const int nc, int **cr);

void add_forces(const struct model m, double **occ, double *forces,
    const int nc, int **cr);

double fermi(const double x);

double dirac(const double x);

double fermi_level(const int ne, double n, const double *e, const double kt,
    double mu);

double grand_potential(const int ne, const double *e, const double kt,
    const double mu);

void occupations(const int ne, const int nspin, double **occ,
    const double *e, double **psi, const double kt, const double mu);
