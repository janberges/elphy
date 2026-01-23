#include "elphy.h"

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
    const char jobz = 'V', uplo = 'U';
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
