#include "elphy.h"

void driver(char *host, double **h, double **h0, double *e, double **occ,
    double **c, double *u, double *forces, const double *forces0,
    double (*tau)[3], const struct model m, const int nc, int **cr,
    const int lwork, double *work) {

    int port, sfd, buf, needinit = 0, havedata = 0;
    char *tmp, header[12];
    const int nat = m.nat * nc;
    int chalen = sizeof(char);
    int intlen = sizeof(int);
    int dbllen = sizeof(double);
    int msglen = 12 * chalen;
    int poslen = 3 * nat * dbllen;
    int matlen = 9 * dbllen;
    double energy, *cell, *positions, *virial;
    int i, j;

    if (!(cell = malloc(matlen)))
        error("No memory for primitive vectors.");
    if (!(positions = malloc(poslen)))
        error("No memory for atomic positions.");
    if (!(virial = calloc(9, sizeof *virial)))
        error("No memory for virial tensor.");

    tmp = strchr(host, ':');

    if (tmp) {
        *tmp = '\0';
        port = atoi(tmp + 1);
    } else
        port = 0;

    if (port)
        sfd = open_inet_socket(host, port);
    else
        sfd = open_unix_socket(host, "/tmp/ipi_");

    for (;;) {
        sread(sfd, header, msglen);

        if (!strncmp(header, "STATUS", 6)) {
            if (needinit)
                swrite(sfd, "NEEDINIT    ", msglen);
            else if (havedata)
                swrite(sfd, "HAVEDATA    ", msglen);
            else
                swrite(sfd, "READY       ", msglen);
        } else if (!strncmp(header, "INIT", 4)) {
            sread(sfd, &buf, intlen); /* replica index */
            sread(sfd, &buf, intlen); /* size of init string */

            if (!(tmp = malloc(buf)))
                error("No memory for init string.");
            sread(sfd, tmp, buf); /* init string */
            free(tmp);

            needinit = 0;
        } else if (!strncmp(header, "POSDATA", 7)) {
            sread(sfd, cell, matlen); /* cell */
            sread(sfd, cell, matlen); /* inverse cell */
            sread(sfd, &buf, intlen); /* number of atoms */
            sread(sfd, positions, poslen); /* positions */

            for (i = 0; i < nat; i++)
                for (j = 0; j < 3; j++)
                    u[3 * i + j] = positions[3 * i + j] - tau[i][j];

            energy = step(h, h0, e, occ, c, u, forces, forces0,
                m, nc, cr, lwork, work);

            havedata = 1;
        } else if (!strncmp(header, "GETFORCE", 8)) {
            swrite(sfd, "FORCEREADY  ", msglen);

            swrite(sfd, &energy, dbllen); /* potential */
            swrite(sfd, &nat, intlen); /* number of atoms */
            swrite(sfd, forces, poslen); /* forces */
            swrite(sfd, virial, matlen); /* virial tensor */

            buf = 1;
            swrite(sfd, &buf, intlen); /* size of extras */
            swrite(sfd, " ", chalen); /* extras */

            havedata = 0;
        } else if (!strncmp(header, "EXIT", 4))
            break;
    }

    free(virial);
    free(positions);
    free(cell);
}
