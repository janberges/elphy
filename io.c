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

int exists(const char *filename) {
    FILE *fp = fopen(filename, "r");

    if (fp) {
        fclose(fp);
        return 1;
    }

    return 0;
}

/* read coupled tight-binding and mass-spring models */

void get_model(const char *filename, struct model *m) {
    FILE *fp;
    int (*r)[3], i, j;
    struct element *t, *k;
    struct vertex *g;
    char *colon;

    fp = fopen(filename, "r");

    if (!fp)
        error("Cannot open %s.", filename);

    if (fscanf(fp, "%s", m->host) != 1)
        error("Invalid i-PI socket address in %s.", filename);

    colon = strchr(m->host, ':');
    if (colon) {
        *colon = '\0';
        m->port = atoi(colon + 1);
    } else
        m->port = 0;

    if (fscanf(fp, "%lf", &m->kt) != 1)
        error("Invalid temperature in %s.", filename);

    if (fscanf(fp, "%lf", &m->umax) != 1)
        error("Invalid radius of random displacements in %s.", filename);

    if (fscanf(fp, "%lf", &m->n) != 1)
        error("Invalid number of electrons per unit cell in %s.", filename);

    if (fscanf(fp, "%d", &m->nel) != 1)
        error("Invalid number of orbitals per unit cell in %s.", filename);

    if (fscanf(fp, "%d", &m->nspin) != 1)
        error("Invalid maximum number of electrons per orbital in %s.",
            filename);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%d", &m->sc[i][j]) != 1)
                error("Invalid supercell vector in %s.", filename);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%lf", &m->uc[i][j]) != 1)
                error("Invalid primitive vector in %s.", filename);

    if (fscanf(fp, "%d", &m->nat) != 1)
        error("Invalid number of atoms per unit cell in %s", filename);
    m->nph = 3 * m->nat;

    m->typ = malloc(m->nat * sizeof *m->typ);
    m->tau = malloc(m->nat * sizeof *m->tau);
    m->fdc = malloc(m->nat * sizeof *m->fdc);

    for (i = 0; i < m->nat; i++) {
        if (fscanf(fp, "%s", m->typ[i]) != 1)
            error("Invalid atom type in %s.", filename);

        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%lf", &m->tau[i][j]) != 1)
                error("Invalid position vector in %s.", filename);

        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%lf", &m->fdc[i][j]) != 1)
                error("Invalid force vector in %s.", filename);
    }

    if (fscanf(fp, "%d", &m->nr) != 1)
        error("Invalid number of lattice vectors in %s.", filename);
    m->r = malloc(m->nr * sizeof *r);

    for (r = m->r; r - m->r < m->nr; r++)
        if (fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2) != 3)
            error("Invalid lattice vector in %s.", filename);

    if (fscanf(fp, "%d", &m->nt) != 1)
        error("Invalid number of hopping parameters in %s.", filename);
    m->t = malloc(m->nt * sizeof *t);

    for (t = m->t; t - m->t < m->nt; t++)
        if (fscanf(fp, "%d %d %d %lf", &t->r, &t->a, &t->b, &t->c) != 4)
            error("Invalid hopping parameter in %s.", filename);

    if (fscanf(fp, "%d", &m->nk) != 1)
        error("Invalid number of interatomic force constants in %s.", filename);
    m->k = malloc(m->nk * sizeof *k);

    for (k = m->k; k - m->k < m->nk; k++)
        if (fscanf(fp, "%d %d %d %lf", &k->r, &k->a, &k->b, &k->c) != 4)
            error("Invalid interatomic force constant in %s.", filename);

    if (fscanf(fp, "%d", &m->ng) != 1)
        error("Invalid number of electron-phonon matrix elements in %s.",
            filename);
    m->g = malloc(m->ng * sizeof *g);

    for (g = m->g; g - m->g < m->ng; g++)
        if (fscanf(fp, "%d %d %d %d %d %lf",
            &g->rph, &g->x, &g->rel, &g->a, &g->b, &g->c) != 6)
                error("Invalid electron-phonon matrix element in %s.",
                    filename);

    fclose(fp);
}

/* get atomic displacements from file with atomic positions in XYZ format */

void get_displ(const char *filename, const int nat, double uc[3][3],
    char (*typ)[3], double (*tau)[3], double *u) {

    const double eps = 1e-5;
    FILE *fp;
    int i, j;
    double r;
    char c[3];

    fp = strcmp(filename, "stdin") ? fopen(filename, "r") : stdin;

    if (!fp)
        error("Cannot open %s.", filename);

    if (fscanf(fp, "%d", &i) != 1)
        error("Invalid number of atoms in %s.", filename);

    if (i != nat)
        error("%d instead of %d atoms in %s.", i, nat, filename);

    if (fscanf(fp, " # CELL{H}:") == EOF)
        error("Unsupported cell definition in %s.", filename);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            if (fscanf(fp, "%lf", &r) != 1)
                error("Invalid cell dimension in %s.", filename);

            if (fabs(r - uc[j][i]) > eps)
                error("Wrong cell dimension in %s.", filename);
        }

    for (i = 0; i < nat; i++) {
        if (fscanf(fp, "%s", c) != 1)
            error("Invalid atom type in %s.", filename);

        if (strcmp(c, typ[i]))
            error("Wrong atom type in %s.", filename);

        for (j = 0; j < 3; j++) {
            if (fscanf(fp, "%lf", &r) != 1)
                error("Invalid position vector %s.", filename);

            u[3 * i + j] = r - tau[i][j];
        }
    }

    if (fp != stdin)
        fclose(fp);
}

/* store positions of displaced atoms in file in XYZ format */

void put_displ(const char *filename, const int nat, double uc[3][3],
    char (*typ)[3], double (*tau)[3], double *u) {

    FILE *fp;
    int i, j;

    fp = strcmp(filename, "stdout") ? fopen(filename, "w") : stdout;

    fprintf(fp, "%d\n", nat);

    fprintf(fp, "# CELL{H}:");

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fprintf(fp, " %.10g", uc[j][i]);

    fprintf(fp, "\n");

    for (i = 0; i < nat; i++) {
        fprintf(fp, "%-3s", typ[i]);
        for (j = 0; j < 3; j++)
            fprintf(fp, " %15.9f", tau[i][j] + u[3 * i + j]);
        fprintf(fp, "\n");
    }

    if (fp != stdout)
        fclose(fp);
}

/* store free energy and forces in file in XYZ format */

void put_force(const char *filename, const int nat, const double energy,
    char (*typ)[3], double (*tau)[3], const double *forces) {

    FILE *fp;
    int i, j;

    fp = strcmp(filename, "stdout") ? fopen(filename, "w") : stdout;

    fprintf(fp, "%d\nfree energy (Ha): %.9f; forces (Ha/bohr):\n", nat, energy);

    for (i = 0; i < nat; i++) {
        fprintf(fp, "%-3s", typ[i]);
        for (j = 0; j < 3; j++)
            fprintf(fp, " %13.9f", forces[3 * i + j]);
        fprintf(fp, "\n");
    }

    if (fp != stdout)
        fclose(fp);
}
