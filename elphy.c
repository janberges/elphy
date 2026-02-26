#include "elphy.h"

#define CC (const char **)
#define CI (const int **)
#define CD (const double **)
#define C3 (const double (*)[3])

int main(const int argc, char **argv) {
    double energy, **h, **h0, *e, **occ, **c, *u, *forces, *forces0, *work, tmp;
    struct model m = {0};
    int nc, nel, nph, nat, **cr, **cells, lwork, info;
    char **typ;
    double (*tau)[3], uc[3][3];

    if (argc > 1 && argc < 5)
        get_model(argv[1], &m);
    else
        error("Usage: elphy <data file> [<socket>|<init file> <radius>]");

    nc = map(m, &cr, &cells);

    nel = m.nel * nc;
    nph = m.nph * nc;
    nat = m.nat * nc;

    if (!(e = malloc(nel * sizeof *e)))
        error("No memory for electron energies.");
    if (!(u = malloc(nph * sizeof *u)))
        error("No memory for atomic displacements.");
    if (!(forces = malloc(nph * sizeof *forces)))
        error("No memory for forces.");
    if (!(forces0 = malloc(nph * sizeof *forces0)))
        error("No memory for force correction.");
    if (!(typ = malloc(nat * sizeof *typ)))
        error("No memory for atom types.");
    if (!(tau = malloc(nat * sizeof *tau)))
        error("No memory for atomic positions.");

    h = matrix(nel);
    h0 = matrix(nel);
    occ = matrix(nel);
    c = matrix(nph);

    lwork = -1;
    dsyev_("V", "U", &nel, *h, &nel, e, &tmp, &lwork, &info);
    lwork = (int) tmp;

    if (!(work = malloc(lwork * sizeof *work)))
        error("No memory for LAPACK work array.");

    populate(h0, m.nel, m.nt, m.t, nc, CI cr);
    populate(c, m.nph, m.nk, m.k, nc, CI cr);

    repeat(uc, typ, tau, (double (*)[3]) forces0, m, nc, CI cells);

    if (argc == 2)
        while (get_displ(nat, CC typ, C3 tau, u) != EOF) {
            energy = step(h, CD h0, e, occ, CD c, u, forces, forces0,
                m, nc, CI cr, lwork, work);

            put_force(nat, CC typ, energy, forces);
        }
    else if (argc == 4) {
        srand(time(NULL));
        random_displacements(nat, u, atof(argv[3]));
        put_displ(argv[2], nat, C3 uc, CC typ, C3 tau, u);

        energy = step(h, CD h0, e, occ, CD c, u, forces, forces0,
            m, nc, CI cr, lwork, work);

        put_force(nat, CC typ, energy, forces);
    } else
        driver(argv[2], h, CD h0, e, occ, CD c, u, forces, forces0, C3 tau,
            m, nc, CI cr, lwork, work);

    free(work);

    free(c);
    free(occ);
    free(h0);
    free(h);

    free(tau);
    free(typ);
    free(forces0);
    free(forces);
    free(u);
    free(e);

    free(cr);
    free(cells);

    free(m.g);
    free(m.k);
    free(m.t);
    free(m.r);
    free(m.fdc);
    free(m.tau);
    free(m.typ);

    return EXIT_SUCCESS;
}

double step(double **h, const double **h0, double *e, double **occ,
    const double **c, const double *u, double *forces, const double *forces0,
    const struct model m, const int nc, const int **cr,
    const int lwork, double *work) {

    double energy;
    static double mu = 0.0;
    const double n = m.n * nc;
    const int nel = m.nel * nc;
    const int nph = m.nph * nc;
    const int inc = 1;
    const double minus = -1.0, zero = 0.0, plus = 1.0;
    int info;

    dsymv_("U", &nph, &minus, *c, &nph, u, &inc, &zero, forces, &inc);

    energy = -0.5 * ddot_(&nph, u, &inc, forces, &inc);

    daxpy_(&nph, &plus, forces0, &inc, forces, &inc);

    memcpy(*h, *h0, nel * nel * sizeof **h);

    perturb(h, u, m, nc, cr);

    dsyev_("V", "U", &nel, *h, &nel, e, work, &lwork, &info);

    mu = fermi_level(n / m.nspin, nel, e, m.kt, mu);

    energy += m.nspin * grand_potential(nel, e, m.kt, mu) + n * mu;

    occupations(nel, e, m.kt, mu, m.nspin, h, occ);

    add_forces(forces, CD occ, m, nc, cr);

    return energy;
}
