/* adapted from i-PI's sockets.c (C) 2013 Joshua More and Michele Ceriotti */
/* connection to INET socket following client example from man getaddrinfo */

#define _POSIX_C_SOURCE 200112L /* beyond ANSI C, see man feature_test_macros */
#include <netdb.h> /* getaddrinfo etc. */
#undef _POSIX_C_SOURCE

#include <sys/socket.h>
#include <sys/un.h> /* UNIX sockets */
#include <unistd.h> /* read and write */
#include <netinet/tcp.h> /* TCP_NODELAY */
#include "elphy.h"

int open_inet_socket(const char *host, int port) {
    const char yes = 1;
    struct addrinfo hints = {0}, *res, *r;
    char service[64];
    int sfd;

    hints.ai_family = AF_UNSPEC; /* IPv4 or IPv6 */
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    hints.ai_protocol = 0; /* any protocol */

    sprintf(service, "%d", port);

    if (getaddrinfo(host, service, &hints, &res))
        error("Could not get address info.");

    for (r = res; r; r = r->ai_next) {
        sfd = socket(r->ai_family, r->ai_socktype, r->ai_protocol);

        if (sfd == -1)
            continue;

        /* see i-PI's sockets.c */
        if (setsockopt(sfd, IPPROTO_TCP, TCP_NODELAY, &yes, sizeof(int)))
            error("Could not set socket option.");

        if (!connect(sfd, r->ai_addr, r->ai_addrlen))
            break;

        close(sfd);
    }

    freeaddrinfo(res);

    if (!r)
        error("Could not connect to INET socket.");

    return sfd;
}

int open_unix_socket(const char *host, const char *prefix) {
    struct sockaddr_un addr = {0};
    int sfd;

    addr.sun_family = AF_UNIX;
    strcpy(addr.sun_path, prefix);
    strcat(addr.sun_path, host);

    sfd = socket(AF_UNIX, SOCK_STREAM, 0);

    if (sfd == -1 || connect(sfd, (struct sockaddr *) &addr, sizeof(addr)))
        error("Could not connect to UNIX socket.");

    return sfd;
}

void swrite(const int sfd, const void *data, const int len) {
    if (write(sfd, data, len) == -1)
        error("Could not write to socket.");
}

void sread(const int sfd, void *data, const int len) {
    int all, new;

    all = 0;
    while (all < len) {
        all += new = read(sfd, (char *) data + all, len - all);
        if (new < 1)
            error("Could not read from socket.");
    }
}
