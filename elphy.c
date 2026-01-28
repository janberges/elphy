#include "elphy.h"

int main(int argc, char **argv) {
    double **h, **c, *e, *u, energy, *forces, *occ, nel, mu = 0.0;
    struct model m;
    int nc, n, nx, i, j, **cr, **cells;
    char (*typ)[3];
    double (*tau)[3], uc[3][3];

    get_model(argc > 1 ? argv[1] : "model.dat", &m);

    nc = map(m, &cr, &cells);

    repeat(m, nc, cells, uc, &typ, &tau);

    nx = m.nph * nc;
    n = m.nel * nc;
    nel = 0.5 * m.n * nc;

    h = matrix(n);
    c = matrix(nx);

    supercell(h, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    u = malloc(nx * sizeof *u);
    get_displ("u.dat", nc * m.nat, uc, typ, tau, u);

    /* put_displ("u_copy.dat", nc * m.nat, uc, typ, tau, u); */

    perturbation(h, m, u, nc, cr);

    e = malloc(n * sizeof *e);
    eigenvalues(n, h, e);

    mu = fermi_level(n, nel, e, m.kt, mu);

    energy = 2.0 * free_energy(n, nel, e, m.kt, mu);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    printf("%.9f\n", energy);

    occ = malloc(n * sizeof *occ);

    for (i = 0; i < n; i++)
        occ[i] = 2.0 * fermi((e[i] - mu) / m.kt);

    forces = malloc(nx * sizeof *forces);
    jacobian(h, m, occ, forces, nc, cr);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            forces[i] += c[i][j] * u[j];

    for (i = 0; i < nx; i++)
        printf("% .9f\n", forces[i]);

    return 0;
}
