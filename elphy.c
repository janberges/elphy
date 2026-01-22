#include <stdio.h>
#include <stdlib.h>

/* declare interface to LAPACK subroutine */

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
    const int *lda, double *w, double *work, const int *lwork, int *info);

double **matrix(const int n);

double *eigenvalues(const int n, double **a);

int main(int argc, char **argv) {
    double **a, *w;
    int n, i, j;

    /* get matrix size from command-line argument */

    n = argc > 1 ? atoi(argv[1]) : 12;

    /* emulate 2D variable-length array (to avoid C99 feature) */

    a = matrix(n);

    /* populate matrix using example of 1D tight-binding Hamiltonian */

    for (i = 0; i < n; i++) {
        j = (i + 1) % n;
        a[i][j] = a[j][i] = -1.0;
    }

    /* diagonalize matrix */

    w = eigenvalues(n, a);

    /* print eigenvalues */

    for (i = 0; i < n; i++)
        printf("% .9f\n", w[i]);

    return 0;
}

double **matrix(const int n) {
    double **a, *data;
    int i;

    a = malloc(n * sizeof data);
    data = calloc(n * n, sizeof *data);

    for (i = 0; i < n; i++)
        a[i] = data + i * n;

    return a;
}

double *eigenvalues(const int n, double **a) {
    const char jobz = 'N', uplo = 'U';
    double *w, *work, lworkopt;
    int lwork, info, i;

    w = malloc(n * sizeof *w);

    for (i = 0; i < 2; i++) {
        if (i) { /* workspace allocation */
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
