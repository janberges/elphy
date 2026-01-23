#include "elphy.h"

int main(int argc, char **argv) {
    const double kt = 0.0019;
    double **h, **c, *e, *u, energy, *forces, *occ, nel = 0.25, mu = 0.0;
    struct model el, ph;
    struct coupling elph;
    int n, nx, nc, i, j;

    get_model("el.dat", &el);
    get_model("ph.dat", &ph);
    get_coupl("elph.dat", &elph);

    /* put_model("el_copy.dat", &el); */

    nc = argc > 1 ? atoi(argv[1]) : 12;
    nx = ph.nb * nc * nc;
    n = el.nb * nc * nc;
    nel *= n;

    h = matrix(n);
    c = matrix(nx);

    supercell(h, el, nc);
    supercell(c, ph, nc);

    u = malloc(nx * sizeof *u);
    get_displ("u.dat", nx, u);

    perturbation(h, elph, u, nc);

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

    forces = jacobian(h, elph, occ, nc);

    for (i = 0; i < nx; i++)
        for (j = 0; j < nx; j++)
            forces[i] += c[i][j] * u[j];

    for (i = 0; i < nx; i++)
        printf("% .9f\n", forces[i]);

    return 0;
}
