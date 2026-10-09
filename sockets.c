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
    int fd;

    hints.ai_family = AF_UNSPEC; /* IPv4 or IPv6 */
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    hints.ai_protocol = 0; /* any protocol */

    if (getaddrinfo(host, port, &hints, &res))
        error("Cannot get address info.");

    for (r = res;; r = r->ai_next ? r->ai_next : (sleep(1), res))
        if ((fd = socket(r->ai_family, r->ai_socktype, r->ai_protocol)) != -1) {
            if (!setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &yes, sizeof yes))
                if (!connect(fd, r->ai_addr, r->ai_addrlen))
                    break;

            close(fd);
        }

    freeaddrinfo(res);

    return fd;
}

int open_unix_socket(const char *host, const char *prefix) {
    struct sockaddr_un addr = {0};
    int fd;

    addr.sun_family = AF_UNIX;
    strncat(addr.sun_path, prefix, sizeof addr.sun_path - 1);
    strncat(addr.sun_path, host, sizeof addr.sun_path - strlen(prefix) - 1);

    if ((fd = socket(AF_UNIX, SOCK_STREAM, 0)) == -1)
        error("Cannot create UNIX socket.");

    while (connect(fd, (struct sockaddr *) &addr, sizeof addr))
        sleep(1);

    return fd;
}

void sread(const int fd, void *data, const int len) {
    int all, new;

    for (all = 0; all < len; all += new)
        if ((new = read(fd, (char *) data + all, len - all)) == -1)
            error("Cannot read from socket.");
}

void swrite(const int fd, const void *data, const int len) {
    if (write(fd, data, len) == -1)
        error("Cannot write to socket.");
}

static void *shmmap(const char *name, const int len) {
    void *addr;
    int fd;

    if ((fd = shm_open(name, O_RDWR, 0)) == -1)
        error("Cannot open shared memory");

    if ((addr = mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0))
            == MAP_FAILED)
        error("Cannot map shared memory");

    close(fd);

    return addr;
}

void *shm_attach(const int fd, const int len) {
    char name[256] = "/";
    int namelen;

    sread(fd, &namelen, sizeof namelen);
    sread(fd, name + 1, namelen);
    name[namelen + 1] = '\0';

    return shmmap(name, len);
}

void shm_detach(void *addr, const int len) {
    munmap(addr, len);
}
