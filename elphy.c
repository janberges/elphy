#include "elphy.h"

#define CC (const char **)
#define CI (const int **)
#define CD (const double **)
#define C3 (const double (*)[3])

int main(const int argc, char **argv) {
    const int inc = 1;
    double energy, energy0, **h, **h0, *e, **occ, **c, *u, *forces, *forces0;
    struct model m = {0};
    int i, n, nc, nel, nph, nat, **cr, **cells, lwork, liwork, *iwork, info;
    char **typ;
    double (*tau)[3], uc[3][3], *work, tmp, *u1, a, b;

    if (argc > 1 && argc < 6)
        get_model(argv[1], &m);
    else
        error("Usage: elphy <data file> "
            "[<socket>|<number> (<radius>|<lower> <upper>)]");

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
    liwork = -1;
    dsyevd_("V", "U", &nel, *h, &nel, e, &tmp, &lwork, &i, &liwork, &info);
    lwork = (int) tmp;
    liwork = i;

    if (!(work = malloc(lwork * sizeof *work)))
        error("No memory for LAPACK work array.");
    if (!(iwork = malloc(liwork * sizeof *iwork)))
        error("No memory for LAPACK integer work array.");

    populate(h0, m.nel, m.nt, m.t, nc, CI cr);
    populate(c, m.nph, m.nk, m.k, nc, CI cr);

    repeat(uc, typ, tau, (double (*)[3]) forces0, m, nc, CI cells);

    if (m.strain) {
        tmp = 1.0 + m.strain;
        info = sizeof uc / sizeof **uc;
        dscal_(&info, &tmp, *uc, &inc);
        dscal_(&nph, &tmp, *tau, &inc);

        strain(h0, m, nc, CI cr);

        energy0 = nc * strain_energy(m);
    } else
        energy0 = 0.0;

    switch (argc) {
    case (2):
        while (get_xyz(nat, CC typ, C3 tau, u) != EOF) {
            energy = step(h, CD h0, e, occ, CD c, u, forces, forces0, energy0,
                m, nc, CI cr, lwork, work, liwork, iwork);

            put_extxyz(nat, C3 uc, CC typ, C3 tau, u, energy, forces);
        }
        break;

    case (3):
        driver(argv[2], h, CD h0, e, occ, CD c, u, forces, forces0, energy0,
            m, nc, CI cr, lwork, work, liwork, iwork, C3 tau);
        break;

    case (4):
        srand(time(NULL));
        n = atoi(argv[2]);

        for (i = 0; i < abs(n); i++) {
            random_displacements(nat, u, atof(argv[3]));

            if (n < 0) {
                put_xyz(nat, C3 uc, CC typ, C3 tau, u, 1);
                continue;
            }

            energy = step(h, CD h0, e, occ, CD c, u, forces, forces0, energy0,
                m, nc, CI cr, lwork, work, liwork, iwork);

            put_extxyz(nat, C3 uc, CC typ, C3 tau, u, energy, forces);
        }
        break;

    case (5):
        u1 = forces; /* use otherwise unused memory */

        for (i = 0; get_xyz(nat, CC typ, C3 tau, u1) != EOF; i++);

        if (!i)
            error("Atomic positions needed.");

        n = atoi(argv[2]);
        a = atof(argv[3]);
        b = atof(argv[4]);

        if (n < 2)
            error("At least two points needed.");

        for (i = 0; i < n; i++) {
            tmp = (a * (n - 1 - i) + b * i) / (n - 1);

            memcpy(u, u1, nph * sizeof *u);
            dscal_(&nph, &tmp, u, &inc);

            put_xyz(nat, C3 uc, CC typ, C3 tau, u, 0);
        }
        break;
    }

    free(iwork);
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
    const double energy0, const struct model m, const int nc, const int **cr,
    const int lwork, double *work, const int liwork, int *iwork) {

    double energy;
    static double mu = 0.0;
    const double n = m.n * nc;
    const int nel = m.nel * nc;
    const int nph = m.nph * nc;
    const int inc = 1;
    const double minus = -1.0, zero = 0.0, plus = 1.0;
    int info;

    dsymv_("U", &nph, &minus, *c, &nph, u, &inc, &zero, forces, &inc);

    energy = energy0 - 0.5 * ddot_(&nph, u, &inc, forces, &inc);

    daxpy_(&nph, &plus, forces0, &inc, forces, &inc);

    memcpy(*h, *h0, nel * nel * sizeof **h);

    perturb(h, u, m, nc, cr);

    dsyevd_("V", "U", &nel, *h, &nel, e, work, &lwork, iwork, &liwork, &info);

    mu = fermi_level(n / m.nspin, nel, e, m.kt, mu);

    energy += m.nspin * grand_potential(nel, e, m.kt, mu) + n * mu;

    occupations(nel, e, m.kt, mu, m.nspin, h, occ);

    add_forces(forces, CD occ, m, nc, cr);

    return energy;
}
