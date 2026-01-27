#include "elphy.h"

int main(int argc, char **argv) {
    const double kt = 0.0019;
    double **h, **c, *e, *u, energy, *forces, *occ, nel = 0.25, mu = 0.0;
    struct model m;
    int n, nx, nc, i, j, **cr;

    get_model("model.dat", &m);

    nc = argc > 1 ? atoi(argv[1]) : 12;
    nx = m.nph * nc * nc;
    n = m.nel * nc * nc;
    nel *= n;

    h = matrix(n);
    c = matrix(nx);

    cr = map(nc, m.nr, m.r);

    supercell(h, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    u = malloc(nx * sizeof *u);
    get_displ("u.dat", nx, u);

    perturbation(h, m, u, nc, cr);

    e = eigenvalues(n, h);

    mu = fermi_level(n, nel, e, kt, mu);

    energy = 2.0 * free_energy(n, nel, e, kt, mu);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            energy += 0.5 * u[i] * c[i][j] * u[j];

    printf("%.9f\n", energy);

    occ = malloc(n * sizeof *occ);

    for (i = 0; i < n; i++)
        occ[i] = 2.0 * fermi((e[i] - mu) / kt);

    forces = jacobian(h, m, occ, nc, cr);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            forces[i] += c[i][j] * u[j];

    for (i = 0; i < nx; i++)
        printf("% .9f\n", forces[i]);

    return 0;
}
