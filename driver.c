#include "elphy.h"

void driver(char *host, double **h, const double **h0, double *e, double **occ,
    const double **c, double *u, double *forces, const double *forces0,
    const double energy0, const struct model m, const int nc, const int **cr,
    const int lwork, double *work, const int liwork, int *iwork,
    const double (*tau)[3]) {

    double energy, *potential = &energy, cell[3][3];
    const double virial[3][3] = {0}, minus = -1.0;
    int sfd, buf, needinit = 0, havedata = 0, shm = 0, attached = 0;
    char *tmp, header[12];
    const int nph = m.nph * nc;
    const int nat = m.nat * nc;
    const int inc = 1;

    tmp = strchr(host, ':');

    if (tmp) {
        *tmp = '\0';
        sfd = open_inet_socket(host, tmp + 1);
    } else {
         tmp = strstr(host, "/shm");

         if (tmp) {
             *tmp = '\0';
             shm = 1;
         }

        sfd = open_unix_socket(host, "/tmp/ipi_");
    }

    for (;;) {
        sread(sfd, header, sizeof header);

        if (!strncmp(header, "STATUS", 6)) {
            if (needinit)
                swrite(sfd, "NEEDINIT    ", sizeof header);
            else if (havedata)
                swrite(sfd, "HAVEDATA    ", sizeof header);
            else
                swrite(sfd, "READY       ", sizeof header);
        } else if (!strncmp(header, "INIT", 4)) {
            sread(sfd, &buf, sizeof buf); /* replica index */
            sread(sfd, &buf, sizeof buf); /* size of init string */

            if (!(tmp = malloc(buf)) && buf)
                error("No memory for init string.");
            sread(sfd, tmp, buf); /* init string */
            free(tmp);

            needinit = 0;
        } else if (!strncmp(header, "POSDATA", 7)) {
            if (!shm) {
                sread(sfd, cell, sizeof cell); /* cell */
                sread(sfd, cell, sizeof cell); /* inverse cell */
            }

            sread(sfd, &buf, sizeof buf); /* number of atoms */

            if (!shm)
                sread(sfd, u, nph * sizeof *u); /* positions */
            else if (!attached) {
                u = shm_attach(sfd, nph * sizeof *u);
                shm_detach(shm_attach(sfd, sizeof cell), sizeof cell);
                shm_detach(shm_attach(sfd, sizeof cell), sizeof cell);
                potential = shm_attach(sfd, sizeof energy);
                forces = shm_attach(sfd, nph * sizeof *forces);
                shm_detach(memset(shm_attach(sfd, sizeof virial), 0,
                    sizeof virial), sizeof virial);

                attached = 1;
            }

            daxpy_(&nph, &minus, *tau, &inc, u, &inc);

            *potential = step(h, h0, e, occ, c, u, forces, forces0, energy0,
                m, nc, cr, lwork, work, liwork, iwork);

            havedata = 1;
        } else if (!strncmp(header, "GETFORCE", 8)) {
            swrite(sfd, "FORCEREADY  ", sizeof header);

            if (!shm) {
                swrite(sfd, potential, sizeof energy);
                swrite(sfd, &nat, sizeof nat);
                swrite(sfd, forces, nph * sizeof *forces);
                swrite(sfd, virial, sizeof virial);
            }

            buf = 1;
            swrite(sfd, &buf, sizeof buf); /* size of extras */
            swrite(sfd, " ", sizeof(char)); /* extras */

            havedata = 0;
        } else if (!strncmp(header, "EXIT", 4)) {
            if (attached) {
                shm_detach(u, nph * sizeof *u);
                shm_detach(potential, sizeof energy);
                shm_detach(forces, nph * sizeof *forces);
            }

            break;
        }
    }
}
