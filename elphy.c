#include <stdio.h>
#include <stdlib.h>

/* declare interface to LAPACK subroutine */

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

double **matrix(const int n);

double *eigenvalues(const int n, double **a);

int main(int argc, char **argv) {
    double **h, *e;
    int n, i, j;

    /* get matrix size from command-line argument */

    n = argc > 1 ? atoi(argv[1]) : 12;

    /* emulate 2D variable-length array (to avoid C99 feature) */

    h = matrix(n);

    /* populate matrix using example of 1D tight-binding Hamiltonian */

    for (i = 0; i < n; i++) {
        j = (i + 1) % n;
        h[i][j] = h[j][i] = -1.0;
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
