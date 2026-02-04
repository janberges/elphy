#include "elphy.h"

void driver(double **h0, double **h, double **c, const struct model m,
    double *u, double *e, double *occ, double *forces, const double *forces0,
    double (*tau)[3], const int nc, int **cr, const int lwork, double *work) {

    const char *sockets_prefix = "/tmp/ipi_";
    int socket, buf, inet = m.port != 0, needinit = 0, havedata = 0;
    char header[12], *initbuffer;
    const int nat = m.nat * nc;
    int chalen = sizeof(char);
    int intlen = sizeof(int);
    int dbllen = sizeof(double);
    int msglen = 12 * chalen;
    int poslen = 3 * nat * dbllen;
    int matlen = 9 * dbllen;
    double energy, *cell, *positions, *virial;
    int i, j;

    cell = malloc(matlen);
    positions = malloc(poslen);
    virial = calloc(9, sizeof *virial);

    open_socket(&socket, &inet, (int *) &m.port, m.host, sockets_prefix);

    for (;;) {
        readbuffer(&socket, header, &msglen);

        if (!strncmp(header, "STATUS", 6)) {
            if (needinit)
                writebuffer(&socket, "NEEDINIT    ", &msglen);
            else if (havedata)
                writebuffer(&socket, "HAVEDATA    ", &msglen);
            else
                writebuffer(&socket, "READY       ", &msglen);
        } else if (!strncmp(header, "INIT", 4)) {
            readbuffer(&socket, &buf, &intlen); /* replica index */
            readbuffer(&socket, &buf, &intlen); /* size of init string */

            initbuffer = malloc(buf);
            readbuffer(&socket, initbuffer, &buf); /* init string */
            free(initbuffer);

            needinit = 0;
        } else if (!strncmp(header, "POSDATA", 7)) {
            readbuffer(&socket, cell, &matlen); /* cell */
            readbuffer(&socket, cell, &matlen); /* inverse cell */
            readbuffer(&socket, &buf, &intlen); /* number of atoms */
            readbuffer(&socket, positions, &poslen); /* positions */

            for (i = 0; i < nat; i++)
                for (j = 0; j < 3; j++)
                    u[3 * i + j] = positions[3 * i + j] - tau[i][j];

            perturbation(h0, h, m, u, nc, cr);

            energy = step(h, c, m, u, e, occ, forces, forces0, nc, cr,
                lwork, work);

            havedata = 1;
        } else if (!strncmp(header, "GETFORCE", 8)) {
            writebuffer(&socket, "FORCEREADY  ", &msglen);

            writebuffer(&socket, &energy, &dbllen); /* potential */
            writebuffer(&socket, &nat, &intlen); /* number of atoms */
            writebuffer(&socket, forces, &poslen); /* forces */
            writebuffer(&socket, virial, &matlen); /* virial tensor */

            buf = 1;
            writebuffer(&socket, &buf, &intlen); /* size of extras */
            writebuffer(&socket, " ", &chalen); /* extras */

            havedata = 0;
        } else if (!strncmp(header, "EXIT", 4))
            break;
    }

    free(virial);
    free(positions);
    free(cell);
}
