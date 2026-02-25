#include "elphy.h"

void error(const char *msg, ...) {
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

    fp = fopen(filename, "r");

    if (!fp)
        error("Cannot open %s.", filename);

    if (fscanf(fp, "%lf", &m->kt) != 1 || m->kt <= 0.0)
        error("Invalid temperature in %s.", filename);

    if (fscanf(fp, "%lf", &m->n) != 1 || m->n < 0.0)
        error("Invalid number of electrons per unit cell in %s.", filename);

    if (fscanf(fp, "%d", &m->nel) != 1 || m->nel < 1)
        error("Invalid number of orbitals per unit cell in %s.", filename);

    if (fscanf(fp, "%d", &m->nspin) != 1 || m->nspin < 1 || m->nspin > 2)
        error("Invalid number of spins per orbital in %s.", filename);

    if (m->n > m->nspin * m->nel)
        error("Electron number too large in %s.", filename);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%d", &m->sc[i][j]) != 1)
                error("Invalid supercell vector in %s.", filename);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%lf", &m->uc[i][j]) != 1)
                error("Invalid primitive vector in %s.", filename);

    if (fscanf(fp, "%d", &m->nat) != 1 || m->nat < 1)
        error("Invalid number of atoms per unit cell in %s", filename);
    m->nph = 3 * m->nat;

    if (!(m->typ = malloc(m->nat * sizeof *m->typ)))
        error("No memory for atom types.");
    if (!(m->tau = malloc(m->nat * sizeof *m->tau)))
        error("No memory for atomic positions.");
    if (!(m->fdc = malloc(m->nat * sizeof *m->fdc)))
        error("No memory for force correction.");

    for (i = 0; i < m->nat; i++) {
        if (fscanf(fp, "%63s", m->typ[i]) != 1)
            error("Invalid atom type in %s.", filename);

        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%lf", &m->tau[i][j]) != 1)
                error("Invalid position vector in %s.", filename);

        for (j = 0; j < 3; j++)
            if (fscanf(fp, "%lf", &m->fdc[i][j]) != 1)
                error("Invalid force vector in %s.", filename);
    }

    if (fscanf(fp, "%d", &m->nr) != 1 || m->nr < 0)
        error("Invalid number of lattice vectors in %s.", filename);
    if (!(m->r = malloc(m->nr * sizeof *r)))
        if (m->nr)
            error("No memory for lattice vectors.");

    for (r = m->r; r - m->r < m->nr; r++)
        if (fscanf(fp, "%d %d %d", *r, *r + 1, *r + 2) != 3)
            error("Invalid lattice vector in %s.", filename);

    if (fscanf(fp, "%d", &m->nt) != 1 || m->nt < 0)
        error("Invalid number of hopping parameters in %s.", filename);
    if (!(m->t = malloc(m->nt * sizeof *t)))
        if (m->nt)
           error("No memory for hopping parameters.");

    for (t = m->t; t - m->t < m->nt; t++) {
        if (fscanf(fp, "%d", &t->r) != 1 || t->r < 0 || t->r >= m->nr)
            error("Invalid lattice-vector index in %s.", filename);
        if (fscanf(fp, "%d", &t->a) != 1 || t->a < 0 || t->a >= m->nel)
            error("Invalid left orbital index in %s.", filename);
        if (fscanf(fp, "%d", &t->b) != 1 || t->b < 0 || t->b >= m->nel)
            error("Invalid right orbital index in %s.", filename);
        if (fscanf(fp, "%lf", &t->c) != 1)
            error("Invalid hopping parameter in %s.", filename);
    }

    if (fscanf(fp, "%d", &m->nk) != 1 || m->nk < 0)
        error("Invalid number of interatomic force constants in %s.", filename);
    if (!(m->k = malloc(m->nk * sizeof *k)))
        if (m->nk)
            error("No memory for interatomic force constants.");

    for (k = m->k; k - m->k < m->nk; k++) {
        if (fscanf(fp, "%d", &k->r) != 1 || k->r < 0 || k->r >= m->nr)
            error("Invalid lattice-vector index in %s.", filename);
        if (fscanf(fp, "%d", &k->a) != 1 || k->a < 0 || k->a >= m->nph)
            error("Invalid left displacement index in %s.", filename);
        if (fscanf(fp, "%d", &k->b) != 1 || k->b < 0 || k->b >= m->nph)
            error("Invalid right displacement index in %s.", filename);
        if (fscanf(fp, "%lf", &k->c) != 1)
            error("Invalid interatomic force constant in %s.", filename);
    }

    if (fscanf(fp, "%d", &m->ng) != 1 || m->ng < 0)
        error("Invalid number of electron-phonon matrix elements in %s.",
            filename);
    if (!(m->g = malloc(m->ng * sizeof *g)))
        if (m->ng)
            error("No memory for electron-phonon matrix elements.");

    for (g = m->g; g - m->g < m->ng; g++) {
        if (fscanf(fp, "%d", &g->rph) != 1 || g->rph < 0 || g->rph >= m->nr)
            error("Invalid lattice-vector index in %s.", filename);
        if (fscanf(fp, "%d", &g->x) != 1 || g->x < 0 || g->x >= m->nph)
            error("Invalid displacement index in %s.", filename);
        if (fscanf(fp, "%d", &g->rel) != 1 || g->rel < 0 || g->rel >= m->nr)
            error("Invalid lattice-vector index in %s.", filename);
        if (fscanf(fp, "%d", &g->a) != 1 || g->a < 0 || g->a >= m->nel)
            error("Invalid left orbital index in %s.", filename);
        if (fscanf(fp, "%d", &g->b) != 1 || g->b < 0 || g->b >= m->nel)
            error("Invalid right orbital index in %s.", filename);
        if (fscanf(fp, "%lf", &g->c) != 1)
            error("Invalid electron-phonon matrix element in %s.", filename);
    }

    fclose(fp);
}

/* input atomic positions in XYZ format and calculate displacements */

int get_displ(const int nat, const char **typ, const double (*tau)[3],
    double *u) {

    int i, j, status;
    double r;
    char c[64];

    if ((status = scanf("%d", &i)) == EOF)
        return EOF;

    if (status != 1)
        error("Invalid number of atoms.");

    if (i != nat)
        error("%d instead of %d atoms.", i, nat);

    if (scanf(" %*[^\n]") != 0)
        error("Invalid comment line.");

    for (i = 0; i < nat; i++) {
        if (scanf("%63s", c) != 1)
            error("Invalid atom type.");

        if (strcmp(c, typ[i]))
            error("Wrong atom type.");

        for (j = 0; j < 3; j++) {
            if (scanf("%lf", &r) != 1)
                error("Invalid position vector.");

            u[3 * i + j] = r - tau[i][j];
        }
    }

    return 0;
}

/* store positions of displaced atoms in file in XYZ format */

void put_displ(const char *filename, const int nat, const double uc[3][3],
    const char **typ, const double (*tau)[3], const double *u) {

    FILE *fp;
    int i, j, width;
    char a[64], *c;

    fp = fopen(filename, "w");

    if (!fp)
        error("Cannot open %s.", filename);

    fprintf(fp, "%d\n", nat);

    fprintf(fp, "# CELL{H}:");

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            sprintf(a, "%.7f", uc[j][i]);
            for (c = a + strlen(a) - 1; c > a && (*c == '0' || *c == '.'); c--)
                *c = '\0';
            fprintf(fp, " %s", a);
        }

    fprintf(fp, "\n");

    width = 0;
    for (i = 0; i < nat; i++)
        if (width < (j = strlen(typ[i])))
            width = j;

    for (i = 0; i < nat; i++) {
        fprintf(fp, "%-*s", width, typ[i]);
        for (j = 0; j < 3; j++)
            fprintf(fp, " %15.9f", tau[i][j] + u[3 * i + j]);
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/* output free energy and forces in XYZ format */

void put_force(const int nat, const char **typ,
    const double energy, const double *forces) {

    int i, j, width;

    printf("%d\nfree energy (Ha): %.9f; forces (Ha/bohr):\n", nat, energy);

    width = 0;
    for (i = 0; i < nat; i++)
        if (width < (j = strlen(typ[i])))
            width = j;

    for (i = 0; i < nat; i++) {
        printf("%-*s", width, typ[i]);
        for (j = 0; j < 3; j++)
            printf(" %13.9f", forces[3 * i + j]);
        printf("\n");
    }
}
