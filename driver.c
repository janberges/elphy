#include "elphy.h"

void driver(char *host, double **h, double **h0, double *e, double **occ,
    double **c, double *u, double *forces, const double *forces0,
    double (*tau)[3], const struct model m, const int nc, int **cr,
    const int lwork, double *work) {

    double energy, cell[3][3];
    const double virial[3][3] = {0};
    int port, sfd, buf, needinit = 0, havedata = 0;
    char *tmp, header[12];
    const int nat = m.nat * nc;
    int i, j;

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

            if (!(tmp = malloc(buf)))
                error("No memory for init string.");
            sread(sfd, tmp, buf); /* init string */
            free(tmp);

            needinit = 0;
        } else if (!strncmp(header, "POSDATA", 7)) {
            sread(sfd, cell, sizeof cell); /* cell */
            sread(sfd, cell, sizeof cell); /* inverse cell */
            sread(sfd, &buf, sizeof buf); /* number of atoms */
            sread(sfd, u, nat * sizeof *tau); /* positions */

            for (i = 0; i < nat; i++)
                for (j = 0; j < 3; j++)
                    u[3 * i + j] -= tau[i][j];

            energy = step(h, h0, e, occ, c, u, forces, forces0,
                m, nc, cr, lwork, work);

            havedata = 1;
        } else if (!strncmp(header, "GETFORCE", 8)) {
            swrite(sfd, "FORCEREADY  ", sizeof header);

            swrite(sfd, &energy, sizeof energy); /* potential */
            swrite(sfd, &nat, sizeof nat); /* number of atoms */
            swrite(sfd, forces, nat * sizeof *tau); /* forces */
            swrite(sfd, virial, sizeof virial); /* virial tensor */

            buf = 1;
            swrite(sfd, &buf, sizeof buf); /* size of extras */
            swrite(sfd, " ", sizeof(char)); /* extras */

            havedata = 0;
        } else if (!strncmp(header, "EXIT", 4))
            break;
    }
}
