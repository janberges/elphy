#include "elphy.h"

int main(int argc, char **argv) {
    double **h0, **h, **c, *e, *u, energy, *forces, *forces0, *cu, *occ;
    struct model m;
    int nc, nel, nph, nat, i, **cr, **cells;
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
    cu = malloc(nph * sizeof *cu);

    h0 = matrix(nel);
    h = matrix(nel);
    c = matrix(nph);

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

            energy = step(h, c, m, u, e, occ, forces, forces0, cu, nc, cr);

            put_force("stdout", nat, energy, typ, tau, forces);
        }
    } else
        driver(h0, h, c, m, u, e, occ, forces, forces0, cu, tau, nc, cr);

    free(*c);
    free(c);
    free(*h);
    free(h);
    free(*h0);
    free(h0);

    free(cu);
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

void random_displacements(const int nat, double *u, double umax) {
    double rho, norm, avg;
    int i, j;

    for (i = 0; i < nat; i++) {
        rho = umax * (double) rand() / (double) RAND_MAX;

        norm = 0.0;

        for (j = 3 * i; j < 3 * i + 3; j++) {
            u[j] = (double) rand();
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
    double *e, double *occ, double *forces, const double *forces0, double *cu,
    const int nc, int **cr) {

    static double mu = 0.0;
    double energy;
    int nel, nph, n, i, j;

    nel = m.nel * nc;
    nph = m.nph * nc;
    n = m.n * nc;

    eigenvalues(nel, h, e);

    mu = fermi_level(nel, m.nspin, n, e, m.kt, mu);

    for (i = 0; i < nph; i++) {
        cu[i] = 0.0;
        for (j = 0; j < nph; j++)
            cu[i] += c[i][j] * u[j];
    }

    energy = free_energy(nel, m.nspin, n, e, m.kt, mu);

    for (i = 0; i < nph; i++)
        energy += 0.5 * u[i] * cu[i];

    occupations(nel, m.nspin, occ, e, m.kt, mu);

    compute_forces(h, m, occ, forces, forces0, nc, cr);

    for (i = 0; i < nph; i++)
        forces[i] -= cu[i];

    return energy;
}
