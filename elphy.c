#include "elphy.h"

int main(int argc, char **argv) {
    double **h, **c, *e, *u, energy, *forces, *occ, nel, mu = 0.0;
    struct model m;
    int n, nx, i, j, **cr;

    get_model(argc > 1 ? argv[1] : "model.dat", &m);

    nx = m.nph * m.nc * m.nc;
    n = m.nel * m.nc * m.nc;
    nel = 0.5 * m.n * m.nc * m.nc;

    h = matrix(n);
    c = matrix(nx);

    cr = map(m.nc, m.nr, m.r);

    supercell(h, m.nel, m.nt, m.t, m.nc, cr);
    supercell(c, m.nph, m.nk, m.k, m.nc, cr);

    u = malloc(nx * sizeof *u);
    get_displ("u.dat", nx, u);

    perturbation(h, m, u, cr);

    e = eigenvalues(n, h);

    mu = fermi_level(n, nel, e, m.kt, mu);

    energy = 2.0 * free_energy(n, nel, e, m.kt, mu);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    printf("%.9f\n", energy);

    occ = malloc(n * sizeof *occ);

    for (i = 0; i < n; i++)
        occ[i] = 2.0 * fermi((e[i] - mu) / m.kt);

    forces = jacobian(h, m, occ, cr);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            forces[i] += c[i][j] * u[j];

    for (i = 0; i < nx; i++)
        printf("% .9f\n", forces[i]);

    return 0;
}
