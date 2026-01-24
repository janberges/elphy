#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

double **matrix(const int n);

int **array_2d(const int rows, const int cols);

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
void put_model(const char *filename, const struct model *m);
void get_coupl(const char *filename, struct coupling *m);
void get_displ(const char *filename, const int nx, double *u);

int **map(const int nc, const int nr, int (*points)[3]);

void supercell(double **h, const struct model m, const int nc, int **cr);

void perturbation(double **h, const struct coupling m, const double *u,
    const int nc, int **cr);

double *jacobian(double **h, const struct coupling m, const double *occ,
    const int nc, int **cr);

double fermi(const double x);
double dirac(const double x);

double fermi_level(const int ne, const double n,
    const double *e, const double kt, double mu);

double free_energy(const int ne, const double n,
    const double *e, const double kt, const double mu);
