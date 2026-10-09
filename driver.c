#include "elphy.h"

void driver(char *host, double **h, const double **h0, double *e, double **occ,
    const double **phi, double *u, double *forces, const double *forces0,
    const double energy0, const struct model m, const int nc, const int **cr,
    const double (*tau)[3]) {

    double energy, *potential = &energy, cell[3][3];
    const double virial[3][3] = {0}, minus = -1.0;
    int fd, buf, needinit = 0, havedata = 0, shm = 0, attached = 0;
    char *tmp, header[12];
    const int nph = m.nph * nc;
    const int nat = m.nat * nc;
    const int inc = 1;

    if (tmp = strchr(host, ':')) {
        *tmp = '\0';
        fd = open_inet_socket(host, tmp + 1);
    } else {
         if (tmp = strstr(host, "/shm")) {
             *tmp = '\0';
             shm = 1;
         }

        fd = open_unix_socket(host, "/tmp/ipi_");
    }

    for (;;) {
        sread(fd, header, sizeof header);

        if (!strncmp(header, "STATUS", 6)) {
            if (needinit)
                swrite(fd, "NEEDINIT    ", sizeof header);
            else if (havedata)
                swrite(fd, "HAVEDATA    ", sizeof header);
            else
                swrite(fd, "READY       ", sizeof header);
        } else if (!strncmp(header, "INIT", 4)) {
            sread(fd, &buf, sizeof buf); /* replica index */
            sread(fd, &buf, sizeof buf); /* size of init string */

            if (!(tmp = malloc(buf)) && buf)
                error("No memory for init string.");
            sread(fd, tmp, buf); /* init string */
            free(tmp);

            needinit = 0;
        } else if (!strncmp(header, "POSDATA", 7)) {
            if (!shm) {
                sread(fd, cell, sizeof cell); /* cell */
                sread(fd, cell, sizeof cell); /* inverse cell */
            }

            sread(fd, &buf, sizeof buf); /* number of atoms */

            if (!shm)
                sread(fd, u, nph * sizeof *u); /* positions */
            else if (!attached) {
                u = shm_attach(fd, nph * sizeof *u);
                shm_detach(shm_attach(fd, sizeof cell), sizeof cell);
                shm_detach(shm_attach(fd, sizeof cell), sizeof cell);
                potential = shm_attach(fd, sizeof energy);
                forces = shm_attach(fd, nph * sizeof *forces);
                shm_detach(memset(shm_attach(fd, sizeof virial), 0,
                    sizeof virial), sizeof virial);

                attached = 1;
            }

            daxpy_(&nph, &minus, *tau, &inc, u, &inc);

            *potential = step(h, h0, e, occ, phi, u, forces, forces0, energy0,
                m, nc, cr);

            havedata = 1;
        } else if (!strncmp(header, "GETFORCE", 8)) {
            swrite(fd, "FORCEREADY  ", sizeof header);

            if (!shm) {
                swrite(fd, potential, sizeof energy);
                swrite(fd, &nat, sizeof nat);
                swrite(fd, forces, nph * sizeof *forces);
                swrite(fd, virial, sizeof virial);
            }

            buf = 1;
            swrite(fd, &buf, sizeof buf); /* size of extras */
            swrite(fd, " ", sizeof(char)); /* extras */

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
