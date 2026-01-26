#include "elphy.h"

/* read tight-binding or mass-spring model */

void get_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3];
    struct element *t;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run data.py first.\n", filename);
        exit(1);
    }

    fscanf(fp, "%d %d %d", &m->nb, &m->nr, &m->nt);

    m->r = malloc(m->nr * sizeof *r);
    m->t = malloc(m->nt * sizeof *t);

    for (r = m->r; r - m->r < m->nr; r++)
        fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2);

    for (t = m->t; t - m->t < m->nt; t++)
        fscanf(fp, "%d %d %d %lf", &t->r, &t->a, &t->b, &t->c);

    fclose(fp);
}

/* read coupling model */

void get_coupl(const char *filename, struct coupling *m) {
    FILE *fp;
    int (*r)[3];
    struct vertex *g;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s. Run data.py first.\n", filename);
        exit(1);
    }

    fscanf(fp, "%d %d %d %d", &m->nph, &m->nel, &m->nr, &m->ng);

    m->r = malloc(m->nr * sizeof *r);
    m->g = malloc(m->ng * sizeof *g);

    for (r = m->r; r - m->r < m->nr; r++)
        fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2);

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
