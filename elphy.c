#include "elphy.h"

int main(int argc, char **argv) {
    double **h, **c, *e, *u, energy, *forces, *occ, n, mu = 0.0;
    struct model m;
    int nc, nel, nph, i, j, **cr, **cells;
    char (*typ)[3];
    double (*tau)[3], uc[3][3];

    get_model(argc > 1 ? argv[1] : "model.dat", &m);

    nc = map(m, &cr, &cells);

    typ = malloc(m.nat * nc * sizeof *typ);
    tau = malloc(m.nat * nc * sizeof *tau);

    repeat(m, nc, cells, uc, typ, tau);

    nel = m.nel * nc;
    nph = m.nph * nc;
    n = 0.5 * m.n * nc;

    h = matrix(nel);
    c = matrix(nph);

    supercell(h, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    u = malloc(nph * sizeof *u);
    get_displ("u.dat", nc * m.nat, uc, typ, tau, u);

    /* put_displ("u_copy.dat", nc * m.nat, uc, typ, tau, u); */

    perturbation(h, m, u, nc, cr);

    e = malloc(nel * sizeof *e);
    eigenvalues(nel, h, e);

    mu = fermi_level(nel, n, e, m.kt, mu);

    energy = 2.0 * free_energy(nel, n, e, m.kt, mu);

    for (i = 0; i < nph; i++)
        for (j = 0; j < nph; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    printf("%.9f\n", energy);

    occ = malloc(nel * sizeof *occ);
    occupations(nel, occ, e, m.kt, mu);

    forces = malloc(nph * sizeof *forces);
    jacobian(h, m, occ, forces, nc, cr);

    for (i = 0; i < nph; i++)
        for (j = 0; j < nph; j++)
            forces[i] += c[i][j] * u[j];

    for (i = 0; i < nph; i++)
        printf("% .9f\n", forces[i]);

    return 0;
}
