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

double *jacobian(double **h, struct coupling m, double *occ, int nc);

double fermi(double x);
double dirac(double x);

double fermi_level(const int ne, const double n,
    const double *e, const double kt, double mu);

double free_energy(const int ne, const double n,
    const double *e, const double kt, const double mu);
