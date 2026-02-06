#include "elphy.h"

int main(int argc, char **argv) {
    double **h0, **h, **c, *e, *u, energy, *forces, *forces0, *occ, tmp, *work;
    struct model m = {0};
    int nc, nel, nph, nat, **cr, **cells, lwork, info;
    char (*typ)[3];
    double (*tau)[3], uc[3][3];

    if (argc > 1 && argc < 5)
        get_model(argv[1], &m);
    else
        error("Usage: elphy <data file> [<socket>|<init file> <radius>]");

    nc = map(m, &cr, &cells);

    nel = m.nel * nc;
    nph = m.nph * nc;
    nat = m.nat * nc;

    if (!(typ = malloc(nat * sizeof *typ)))
        error("No memory for atom types.");
    if (!(tau = malloc(nat * sizeof *tau)))
        error("No memory for atomic positions.");
    if (!(u = malloc(nph * sizeof *u)))
        error("No memory for atomic displacements.");
    if (!(e = malloc(nel * sizeof *e)))
        error("No memory for electron energies.");
    if (!(occ = malloc(nel * sizeof *occ)))
        error("No memory for electron occupations.");
    if (!(forces = malloc(nph * sizeof *forces)))
        error("No memory for forces.");
    if (!(forces0 = malloc(nph * sizeof *forces0)))
        error("No memory for force correction.");

    h0 = matrix(nel);
    h = matrix(nel);
    c = matrix(nph);

    lwork = -1;
    dsyev_("V", "U", &nel, *h, &nel, e, &tmp, &lwork, &info);
    lwork = (int) tmp;

    if (!(work = malloc(lwork * sizeof *work)))
        error("No memory for LAPACK work array.");

    supercell(h0, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    repeat(m, nc, cells, uc, typ, tau, (double (*)[3]) forces0);

    if (argc == 2)
        while (get_displ("stdin", nat, uc, typ, tau, u) != EOF) {
            energy = step(h0, h, c, m, u, e, occ, forces, forces0, nc, cr,
                lwork, work);

            put_force("stdout", nat, energy, typ, tau, forces);
        }
    else if (argc == 4) {
        srand(time(NULL));
        random_displacements(nat, u, atof(argv[3]));
        put_displ(argv[2], nat, uc, typ, tau, u);

        energy = step(h0, h, c, m, u, e, occ, forces, forces0, nc, cr,
            lwork, work);

        put_force("stdout", nat, energy, typ, tau, forces);
    } else
        driver(h0, h, c, m, u, e, occ, forces, forces0, tau, nc, cr,
            lwork, work, argv[2]);

    free(work);

    free(*c);
    free(c);
    free(*h);
    free(h);
    free(*h0);
    free(h0);

    free(forces0);
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
    free(m.fdc);
    free(m.tau);
    free(m.typ);

    return 0;
}

double step(double **h0, double **h, double **c, const struct model m,
    const double *u, double *e, double *occ, double *forces,
    const double *forces0, const int nc, int **cr, const int lwork,
    double *work) {

    const int inc = 1;
    const double minus = -1.0, zero = 0.0, plus = 1.0;
    static double mu = 0.0;
    double energy;
    const double n = m.n * nc;
    const int nel = m.nel * nc;
    const int nph = m.nph * nc;
    int info;

    perturbation(h0, h, m, u, nc, cr);

    dsyev_("V", "U", &nel, *h, &nel, e, work, &lwork, &info);

    mu = fermi_level(nel, n / m.nspin, e, m.kt, mu);

    dsymv_("U", &nph, &minus, *c, &nph, u, &inc, &zero, forces, &inc);

    energy = -0.5 * ddot_(&nph, u, &inc, forces, &inc);

    energy += m.nspin * grand_potential(nel, e, m.kt, mu) + n * mu;

    occupations(nel, m.nspin, occ, e, m.kt, mu);

    add_forces(h, m, occ, forces, nc, cr);

    daxpy_(&nph, &plus, forces0, &inc, forces, &inc);

    return energy;
}
