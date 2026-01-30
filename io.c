#include "elphy.h"

void error(char *msg, ...) {
    va_list args;
    char fmt[256] = "elphy error: ";

    strcat(fmt, msg);
    strcat(fmt, "\n");

    va_start(args, msg);
    vfprintf(stderr, fmt, args);
    va_end(args);

    exit(EXIT_FAILURE);
}

/* read coupled tight-binding and mass-spring models */

void get_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3], i, j;
    struct element *t, *k;
    struct vertex *g;
    char *colon;

    fp = fopen(filename, "r");

    if (fp == NULL)
        error("Cannot open %s. Run input.py first.", filename);

    fscanf(fp, "%s", m->host);
    colon = strchr(m->host, ':');
    if (colon) {
        *colon = '\0';
        m->port = atoi(colon + 1);
    } else
        m->port = 0;

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

    m->typ = malloc(m->nat * sizeof *m->typ);
    m->tau = malloc(m->nat * sizeof *m->tau);

    for (i = 0; i < m->nat; i++) {
        fscanf(fp, "%s", m->typ[i]);
        for (j = 0; j < 3; j++)
            fscanf(fp, "%lf", &m->tau[i][j]);
    }

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

/* get atomic displacements from file with atomic positions in XYZ format */

void get_displ(const char *filename, const int nat,
    double uc[3][3], char (*typ)[3], double (*tau)[3], double *u) {

    const double eps = 1e-5;
    FILE *fp;
    int i, j;
    double r;
    char c[3];

    fp = fopen(filename, "r");

    if (fp == NULL)
        error("Cannot open %s. Run test.py first.", filename);

    fscanf(fp, "%d", &i);

    if (i != nat)
        error("%d instead of %d atoms in %s.", i, nat, filename);

    if (fscanf(fp, " # CELL{H}:") == EOF)
        error("Unsupported cell definition in %s.", filename);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            fscanf(fp, "%lf", &r);
            if (fabs(r - uc[j][i]) > eps)
                error("Wrong cell dimension in %s.", filename);
        }

    for (i = 0; i < nat; i++) {
        fscanf(fp, "%s", c);

        if (strcmp(c, typ[i]))
            error("Wrong atom type in %s.", filename);

        for (j = 0; j < 3; j++) {
            fscanf(fp, "%lf", &r);
            u[3 * i + j] = r - tau[i][j];
        }
    }

    fclose(fp);
}

/* store positions of displaced atoms in file in XYZ format */

void put_displ(const char *filename, const int nat,
    double uc[3][3], char (*typ)[3], double (*tau)[3], double *u) {

    FILE *fp;
    int i, j;

    fp = fopen(filename, "w");

    fprintf(fp, "%d\n", nat);

    fprintf(fp, "# CELL{H}:");

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fprintf(fp, " %.10g", uc[j][i]);

    fprintf(fp, "\n");

    for (i = 0; i < nat; i++) {
        fprintf(fp, "%8s", typ[i]);
        for (j = 0; j < 3; j++)
            fprintf(fp, " %15.9f", tau[i][j] + u[3 * i + j]);
        fprintf(fp, "\n");
    }

    fclose(fp);
}
