#include "elphy.h"

/* emulate 2D variable-length array (to avoid C99 feature) */

double **matrix(const int n) {
    double **a, *data;
    int row;

    if (!(a = malloc(n * (sizeof data + n * sizeof *data))))
        if (n)
            error("No memory for %d x %d doubles.", n, n);

    data = (double *) (a + n);

    memset(data, 0, n * n * sizeof *data);

    for (row = 0; row < n; row++)
        a[row] = data + row * n;

    return a;
}

/* alternative version for non-square integer arrays */

int **array_2d(const int rows, const int cols) {
    int **a, *data;
    int row;

    if (!(a = malloc(rows * (sizeof data + cols * sizeof *data))))
        if (rows)
            error("No memory for %d x %d integers.", rows, cols);

    data = (int *) (a + rows);

    memset(data, 0, rows * cols * sizeof *data);

    for (row = 0; row < rows; row++)
        a[row] = data + row * cols;

    return a;
}
