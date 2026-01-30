/* adapted from i-PI's sockets.c (C) 2013 Joshua More and Michele Ceriotti */
/* connection to INET socket following client example from man getaddrinfo */

#define _POSIX_C_SOURCE 200112L /* beyond ANSI C, see man feature_test_macros */
#include <netdb.h> /* getaddrinfo etc. */
#include <sys/socket.h>
#include <sys/un.h> /* UNIX sockets */
#include <unistd.h> /* read and write */
#include <netinet/tcp.h> /* TCP_NODELAY */
#include "elphy.h"

void open_socket(int *psockfd, int *inet, int *port, const char *host,
    const char *sockets_prefix) {

    const char yes = 1;
    struct addrinfo hints = {0}, *res, *r;
    struct sockaddr_un addr = {0};
    char service[64];
    int sfd;

    if (*inet) {
        hints.ai_family = AF_UNSPEC; /* IPv4 or IPv6 */
        hints.ai_socktype = SOCK_STREAM;
        hints.ai_flags = AI_PASSIVE;
        hints.ai_protocol = 0; /* any protocol */

        sprintf(service, "%d", *port);

        if (getaddrinfo(host, service, &hints, &res)) {
            fprintf(stderr, "Could not get address info.\n");
            exit(1);
        }

        for (r = res; r; r = r->ai_next) {
            sfd = socket(r->ai_family, r->ai_socktype, r->ai_protocol);

            if (sfd == -1)
                continue;

            /* see i-PI's sockets.c */
            if (setsockopt(sfd, IPPROTO_TCP, TCP_NODELAY, &yes, sizeof(int))) {
                fprintf(stderr, "Could not set socket option.\n");
                exit(1);
            }

            if (!connect(sfd, r->ai_addr, r->ai_addrlen))
                break;

            close(sfd);
        }

        freeaddrinfo(res);

        if (!r) {
            fprintf(stderr, "Could not connect to INET socket.\n");
            exit(1);
        }
    } else {
        addr.sun_family = AF_UNIX;
        strcpy(addr.sun_path, sockets_prefix);
        strcat(addr.sun_path, host);

        sfd = socket(AF_UNIX, SOCK_STREAM, 0);

        if (sfd == -1 || connect(sfd, (struct sockaddr *) &addr, sizeof(addr))) {
            fprintf(stderr, "Could not connect to UNIX socket.\n");
            exit(1);
        }
    }

    *psockfd = sfd;
}

void writebuffer(int *psockfd, const void *data, int *plen) {
    if (write(*psockfd, data, *plen) == -1) {
        fprintf(stderr, "Could not write to socket.\n");
        exit(1);
    }
}

void readbuffer(int *psockfd, void *data, int *plen) {
    int len, new;

    len = 0;
    while (len < *plen) {
        len += new = read(*psockfd, (char *) data + len, *plen - len);
        if (new < 1) {
            fprintf(stderr, "Could not read from socket.\n");
            exit(1);
        }
    }
}
