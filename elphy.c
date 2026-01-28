#include "elphy.h"

int main(int argc, char **argv) {
    double **h, **c, *e, *u, energy, *forces, *occ, n, mu = 0.0;
    struct model m;
    int nc, nel, nph, nat, i, j, **cr, **cells;
    char (*typ)[3];
    double (*tau)[3], uc[3][3];

    get_model(argc > 1 ? argv[1] : "model.dat", &m);

    nc = map(m, &cr, &cells);

    nel = m.nel * nc;
    nph = m.nph * nc;
    nat = m.nat * nc;
    n = 0.5 * m.n * nc;

    typ = malloc(nat * sizeof *typ);
    tau = malloc(nat * sizeof *tau);
    u = malloc(nph * sizeof *u);
    e = malloc(nel * sizeof *e);
    occ = malloc(nel * sizeof *occ);
    forces = malloc(nph * sizeof *forces);

    h = matrix(nel);
    c = matrix(nph);

    supercell(h, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    repeat(m, nc, cells, uc, typ, tau);

    get_displ("u.dat", nat, uc, typ, tau, u);

    /* put_displ("u_copy.dat", nat, uc, typ, tau, u); */

    perturbation(h, m, u, nc, cr);

    eigenvalues(nel, h, e);

    mu = fermi_level(nel, n, e, m.kt, mu);

    energy = 2.0 * free_energy(nel, n, e, m.kt, mu);

    for (i = 0; i < nph; i++)
        for (j = 0; j < nph; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    printf("%.9f\n", energy);

    occupations(nel, occ, e, m.kt, mu);

    jacobian(h, m, occ, forces, nc, cr);

    for (i = 0; i < nph; i++)
        for (j = 0; j < nph; j++)
            forces[i] += c[i][j] * u[j];

    for (i = 0; i < nph; i++)
        printf("% .9f\n", forces[i]);

    free(*c);
    free(c);
    free(*h);
    free(h);

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
