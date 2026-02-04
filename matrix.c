#include "elphy.h"

/* emulate 2D variable-length array (to avoid C99 feature) */

double **matrix(const int n) {
    double **a, *data;
    int row;

    a = malloc(n * sizeof data);
    data = calloc(n * n, sizeof *data);

    for (row = 0; row < n; row++)
        a[row] = data + row * n;

    return a; /* free data using free(*a) before free(a) */
}

/* alternative version for non-square integer arrays */

int **array_2d(const int rows, const int cols) {
    int **a, *data;
    int row;

    a = malloc(rows * sizeof data);
    data = calloc(rows * cols, sizeof *data);

    for (row = 0; row < rows; row++)
        a[row] = data + row * cols;

    return a; /* free data using free(*a) before free(a) */
}
