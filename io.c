#include "elphy.h"

/* read coupled tight-binding and mass-spring models */

void get_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3], i, j;
    struct element *t, *k;
    struct vertex *g;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run data.py first.\n", filename);
        exit(1);
    }

    fscanf(fp, "%lf", &m->kt);
    fscanf(fp, "%lf", &m->n);
    fscanf(fp, "%d", &m->nel);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fscanf(fp, "%d", &m->sc[i][j]);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fscanf(fp, "%lf", &m->uc[i][j]);

    fscanf(fp, "%d", &m->nat);
    m->nph = 3 * m->nat;

    m->tau = malloc(m->nat * sizeof *m->tau);

    for (i = 0; i < m->nat; i++)
        for (j = 0; j < 3; j++)
            fscanf(fp, "%lf", &m->tau[i][j]);

    fscanf(fp, "%d", &m->nr);
    m->r = malloc(m->nr * sizeof *r);

    for (r = m->r; r - m->r < m->nr; r++)
        fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2);

    fscanf(fp, "%d", &m->nt);
    m->t = malloc(m->nt * sizeof *t);

    for (t = m->t; t - m->t < m->nt; t++)
        fscanf(fp, "%d %d %d %lf", &t->r, &t->a, &t->b, &t->c);

    fscanf(fp, "%d", &m->nk);
    m->k = malloc(m->nk * sizeof *k);

    for (k = m->k; k - m->k < m->nk; k++)
        fscanf(fp, "%d %d %d %lf", &k->r, &k->a, &k->b, &k->c);

    fscanf(fp, "%d", &m->ng);
    m->g = malloc(m->ng * sizeof *g);

    for (g = m->g; g - m->g < m->ng; g++)
        fscanf(fp, "%d %d %d %d %d %lf",
            &g->rph, &g->x, &g->rel, &g->a, &g->b, &g->c);

    fclose(fp);
}

/* read atomic displacements */

void get_displ(const char *filename, const int nx, double *u) {
    FILE *fp;
    int i;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run test.py first.\n", filename);
        exit(1);
    }

    for (i = 0; i < nx; i++)
        fscanf(fp, "%lf", u++);

    fclose(fp);
}
