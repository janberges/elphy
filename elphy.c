#include "elphy.h"

int main(int argc, char **argv) {
    double **h0, **h, **c, *e, *u, energy, *forces, *forces0, *occ;
    struct model m;
    int nc, nel, nph, nat, **cr, **cells;
    char (*typ)[3];
    double (*tau)[3], uc[3][3];

    get_model(argc > 1 ? argv[1] : "input.dat", &m);

    nc = map(m, &cr, &cells);

    nel = m.nel * nc;
    nph = m.nph * nc;
    nat = m.nat * nc;

    typ = malloc(nat * sizeof *typ);
    tau = malloc(nat * sizeof *tau);
    u = malloc(nph * sizeof *u);
    e = malloc(nel * sizeof *e);
    occ = malloc(nel * sizeof *occ);
    forces = malloc(nph * sizeof *forces);
    forces0 = calloc(nph, sizeof *forces);

    h0 = matrix(nel);
    h = matrix(nel);
    c = matrix(nph);

    supercell(h0, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    memcpy(*h, *h0, nel * nel * sizeof **h);

    energy = step(h, c, m, u, e, occ, forces, forces0, nc, cr);

    repeat(m, nc, cells, uc, typ, tau);

    if (!strcmp(m.host, "none")) {
        get_displ("input.xyz", nat, uc, typ, tau, u);

        perturbation(h0, h, m, u, nc, cr);

        energy = step(h, c, m, u, e, occ, forces, forces0, nc, cr);

        put_force("stdout", nat, energy, typ, tau, forces);
    } else
        driver(h0, h, c, m, u, e, occ, forces, forces0, tau, nc, cr);

    free(*c);
    free(c);
    free(*h);
    free(h);
    free(*h0);
    free(h0);

    free(forces);
    free(occ);
    free(e);
    free(u);
    free(tau);
    free(typ);

    free(*cr);
    free(cr);
    free(*cells);
    free(cells);

    free(m.g);
    free(m.k);
    free(m.t);
    free(m.r);
    free(m.tau);
    free(m.typ);

    return 0;
}

double step(double **h, double **c, const struct model m, const double *u,
    double *e, double *occ, double *forces, const double *forces0,
    const int nc, int **cr) {

    static double mu = 0.0;
    double energy;
    int nel, nph, n, i, j;

    nel = m.nel * nc;
    nph = m.nph * nc;
    n = 0.5 * m.n * nc;

    eigenvalues(nel, h, e);

    mu = fermi_level(nel, n, e, m.kt, mu);

    energy = 2.0 * free_energy(nel, n, e, m.kt, mu);

    for (i = 0; i < nph; i++)
        for (j = 0; j < nph; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    occupations(nel, occ, e, m.kt, mu);

    jacobian(h, m, occ, forces, forces0, nc, cr);

    for (i = 0; i < nph; i++) {
        for (j = 0; j < nph; j++)
            forces[i] += c[i][j] * u[j];
        forces[i] *= -1;
    }

    return energy;
}
