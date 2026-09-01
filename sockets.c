/* adapted from i-PI's sockets.c (C) 2013 Joshua More and Michele Ceriotti */
/* connection to INET socket following client example from man getaddrinfo */

#define _POSIX_C_SOURCE 200112L /* beyond ANSI C, see man feature_test_macros */
#include <netdb.h> /* getaddrinfo etc. */
#undef _POSIX_C_SOURCE

#include <sys/socket.h>
#include <sys/un.h> /* UNIX sockets */
#include <sys/mman.h> /* memory management */
#include <unistd.h> /* read and write */
#include <netinet/tcp.h> /* TCP_NODELAY */
#include <fcntl.h> /* O_RDWR */
#include "elphy.h"

int open_inet_socket(const char *host, const char *port) {
    const int yes = 1;
    struct addrinfo hints = {0}, *res, *r;
    int sfd;

    hints.ai_family = AF_UNSPEC; /* IPv4 or IPv6 */
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    hints.ai_protocol = 0; /* any protocol */

    if (getaddrinfo(host, port, &hints, &res))
        error("Cannot get address info.");

    for (r = res; r; r = r->ai_next) {
        sfd = socket(r->ai_family, r->ai_socktype, r->ai_protocol);

        if (sfd == -1)
            continue;

        /* see i-PI's sockets.c */
        if (setsockopt(sfd, IPPROTO_TCP, TCP_NODELAY, &yes, sizeof yes))
            error("Cannot set socket option.");

        if (!connect(sfd, r->ai_addr, r->ai_addrlen))
            break;

        close(sfd);
    }

    freeaddrinfo(res);

    if (!r)
        error("Cannot connect to INET socket.");

    return sfd;
}

int open_unix_socket(const char *host, const char *prefix) {
    struct sockaddr_un addr = {0};
    int sfd;

    addr.sun_family = AF_UNIX;
    strncat(addr.sun_path, prefix, sizeof addr.sun_path - 1);
    strncat(addr.sun_path, host, sizeof addr.sun_path - strlen(prefix) - 1);

    sfd = socket(AF_UNIX, SOCK_STREAM, 0);

    if (sfd == -1 || connect(sfd, (struct sockaddr *) &addr, sizeof addr))
        error("Cannot connect to UNIX socket.");

    return sfd;
}

void sread(const int sfd, void *data, const int len) {
    int all, new;

    all = 0;
    while (all < len) {
        all += new = read(sfd, (char *) data + all, len - all);
        if (new < 1)
            error("Cannot read from socket.");
    }
}

void swrite(const int sfd, const void *data, const int len) {
    if (write(sfd, data, len) == -1)
        error("Cannot write to socket.");
}

static void *shmmap(const char *name, const int len) {
    void *addr;
    int mfd;

    mfd = shm_open(name, O_RDWR, 0666);

    if (mfd == -1)
        error("Cannot open shared memory");

    addr = mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_SHARED, mfd, 0);

    if (addr == MAP_FAILED)
        error("Cannot map shared memory");

    close(mfd);

    return addr;
}

void *shm_attach(const int sfd, const int len) {
    char name[256] = "/";
    int namelen;

    sread(sfd, &namelen, sizeof namelen);
    sread(sfd, name + 1, namelen);
    name[namelen + 1] = '\0';

    return shmmap(name, len);
}

void shm_detach(void *addr, const int len) {
    munmap(addr, len);
}
