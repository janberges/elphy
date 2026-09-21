#include "elphy.h"

#define CC (const char **)
#define CI (const int **)
#define CD (const double **)
#define C3 (const double (*)[3])

static int lwork, liwork, *iwork;
static double *work;

int main(const int argc, char **argv) {
    const int inc = 1;
    double **h, **h0, *e, **occ, **c, *u, *forces, *forces0, energy, energy0,
        (*tau)[3], uc[3][3], tmp, *u1, a, b, dt, damp, s, *swap;
    struct model m = {0};
    int i, j, n, nc, nel, nph, nat, **cr, **cells, info, stride;
    char **typ, *match;

    if (argc > 1 && argc < 7)
        get_model(argv[1], &m);
    else
        error("Usage: elphy <data file> [<socket>|<number> "
            "(<radius>|<lower> <upper>|<dt> <damp> <vmax>)]");

    nc = map(m, &cr, &cells);

    nel = m.nel * nc;
    nph = m.nph * nc;
    nat = m.nat * nc;

    if (!(e = malloc(nel * sizeof *e)))
        error("No memory for electron energies.");
    if (!(u = malloc(nph * sizeof *u)))
        error("No memory for atomic displacements.");
    if (!(u1 = malloc(nph * sizeof *u1)))
        error("No memory for original displacements.");
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

    if (m.strain) {
        strain(h0, m, nc, CI cr);

        energy0 = nc * strain_model(&m);
    } else
        energy0 = 0.0;

    repeat(uc, typ, tau, (double (*)[3]) forces0, m, nc, CI cells);

    switch (argc) {
    case (2):
        while (get_xyz(nat, CC typ, C3 tau, u) != EOF) {
            energy = step(h, CD h0, e, occ, CD c, u, forces, forces0, energy0,
                m, nc, CI cr);

            put_extxyz(nat, C3 uc, CC typ, C3 tau, u, energy, forces);
        }
        break;

    case (3):
        driver(argv[2], h, CD h0, e, occ, CD c, u, forces, forces0, energy0,
            m, nc, CI cr, C3 tau);
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
                m, nc, CI cr);

            put_extxyz(nat, C3 uc, CC typ, C3 tau, u, energy, forces);
        }
        break;

    case (5):
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

    case (6):
        match = strchr(argv[2], ':');

        if (match) {
            *match = '\0';
            if ((stride = atoi(match + 1)) < 1)
                error("Stride must be at least one.");
        } else
            stride = 1;

        n = atoi(argv[2]);

        if (!(dt = atof(argv[3])))
            error("Time step must be nonzero.");

        damp = 0.5 * atof(argv[4]) * dt;

        if (damp < 0.0) {
            damp *= -1.0;
            s = 0.0;
        } else
            s = 2.0 / dt * sqrt(damp * m.kt);

        memset(u1, 0, nph * sizeof *u1);
        random_displacements(nat, u, atof(argv[5]) * dt);

        a = 2.0 / (1.0 + damp);
        b = (damp - 1.0) / (1.0 + damp);
        tmp = dt * dt / (1.0 + damp);

        for (i = 0; i < n; i++) {
            energy = step(h, CD h0, e, occ, CD c, u, forces, forces0, energy0,
                m, nc, CI cr);

            if (!(i % stride))
                put_extxyz(nat, C3 uc, CC typ, C3 tau, u, energy, forces);

            for (j = 0; j < nph; j++) {
                if (s)
                    forces[j] += s * sqrt(m.mass[j / 3 % m.nat]) * box_muller();

                forces[j] /= m.mass[j / 3 % m.nat];
            }

            dscal_(&nph, &b, u1, &inc);
            daxpy_(&nph, &a, u, &inc, u1, &inc);
            daxpy_(&nph, &tmp, forces, &inc, u1, &inc);

            fixcom(nat, u1);

            swap = u;
            u = u1;
            u1 = swap;
        }
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
    free(u1);
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
    free(m.mass);
    free(m.typ);

    return EXIT_SUCCESS;
}

double step(double **h, const double **h0, double *e, double **occ,
    const double **c, const double *u, double *forces, const double *forces0,
    const double energy0, const struct model m, const int nc, const int **cr) {

    double energy;
    static double mu = 0.0;
    const double n = m.n * nc;
    const int nel = m.nel * nc;
    const int nph = m.nph * nc;
    const int inc = 1;
    const double minus = -1.0, plus = 1.0;
    int info;

    memcpy(forces, forces0, nph * sizeof *forces);

    dsymv_("U", &nph, &minus, *c, &nph, u, &inc, &plus, forces, &inc);

    energy = energy0 - 0.5 * ddot_(&nph, u, &inc, forces, &inc);
    energy = energy - 0.5 * ddot_(&nph, u, &inc, forces0, &inc);

    memcpy(*h, *h0, nel * nel * sizeof **h);

    perturb(h, u, m, nc, cr);

    dsyevd_("V", "U", &nel, *h, &nel, e, work, &lwork, iwork, &liwork, &info);

    mu = fermi_level(n / m.nspin, nel, e, m.kt, mu);

    energy += m.nspin * grand_potential(nel, e, m.kt, mu) + n * mu;

    occupations(nel, e, m.kt, mu, m.nspin, h, occ);

    add_forces(forces, CD occ, m, nc, cr);

    return energy;
}
