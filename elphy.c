#include "elphy.h"

int main(int argc, char **argv) {
    double **h0, **h, **c, *e, *u, energy, *forces, *forces0, *occ, tmp, *work;
    struct model m = {0};
    int nc, nel, nph, nat, i, **cr, **cells, lwork, info;
    char (*typ)[3];
    double (*tau)[3], uc[3][3];

    if (argc > 1)
        get_model(argv[1], &m);
    else
        error("Argument missing. Usage: elphy input.dat [input.xyz ...]");

    nc = map(m, &cr, &cells);

    nel = m.nel * nc;
    nph = m.nph * nc;
    nat = m.nat * nc;

    typ = malloc(nat * sizeof *typ);
    tau = malloc(nat * sizeof *tau);
    u = malloc(nph * sizeof *u);
    e = malloc(nel * sizeof *e);
    occ = malloc(nel * sizeof *occ);
    forces = malloc(nph * sizeof *forces);
    forces0 = malloc(nph * sizeof *forces0);

    h0 = matrix(nel);
    h = matrix(nel);
    c = matrix(nph);

    lwork = -1;
    dsyev_("V", "U", &nel, *h, &nel, e, &tmp, &lwork, &info);
    lwork = (int) tmp;

    work = malloc(lwork * sizeof *work);

    supercell(h0, m.nel, m.nt, m.t, nc, cr);
    supercell(c, m.nph, m.nk, m.k, nc, cr);

    repeat(m, nc, cells, uc, typ, tau, (double (*)[3]) forces0);

    if (argc > 2) {
        for (i = 2; i < argc; i++) {
            if (exists(argv[i]))
                get_displ(argv[i], nat, uc, typ, tau, u);
            else {
                random_displacements(nat, u, m.umax);
                put_displ(argv[i], nat, uc, typ, tau, u);
            }

            perturbation(h0, h, m, u, nc, cr);

            energy = step(h, c, m, u, e, occ, forces, forces0, nc, cr,
                lwork, work);

            put_force("stdout", nat, energy, typ, tau, forces);
        }
    } else
        driver(h0, h, c, m, u, e, occ, forces, forces0, tau, nc, cr,
            lwork, work);

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

/* MINSTD using double instead of long long int (to avoid C99 feature) */

double minstd(void) {
    const double a = 48271.0, m = 2147483647.0;
    static double i = 1.0;

    i = fmod(i * a, m);

    return i / m;
}

void random_displacements(const int nat, double *u, double umax) {
    double rho, norm, avg;
    int i, j;

    for (i = 0; i < nat; i++) {
        rho = umax * minstd();

        norm = 0.0;

        for (j = 3 * i; j < 3 * i + 3; j++) {
            u[j] = 1.0 - 2.0 * minstd();
            norm += u[j] * u[j];
        }

        norm = sqrt(norm);

        if (norm != 0.0)
            for (j = 0; j < 3; j++)
                u[3 * i + j] *= rho / norm;
    }

    for (j = 0; j < 3; j++) {
        avg = 0.0;

        for (i = 0; i < nat; i++)
            avg += u[3 * i + j];

        avg /= nat;

        for (i = 0; i < nat; i++)
            u[3 * i + j] -= avg;
    }
}

double step(double **h, double **c, const struct model m, const double *u,
    double *e, double *occ, double *forces, const double *forces0, const int nc,
    int **cr, const int lwork, double *work) {

    const int inc = 1;
    const double minus = -1.0, zero = 0.0, plus = 1.0;
    static double mu = 0.0;
    double energy;
    const double n = m.n * nc;
    const int nel = m.nel * nc;
    const int nph = m.nph * nc;
    int info;

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
