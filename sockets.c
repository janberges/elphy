/* adapted from i-PI's sockets.c (C) 2013 Joshua More and Michele Ceriotti */
/* connection to INET socket following client example from man getaddrinfo */

#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#include <windows.h> /* Sleep */

#define sleep(seconds) Sleep(1000 * seconds)
#define close(fd) closesocket(fd)

static SOCKET fd;
#else
#define _POSIX_C_SOURCE 200112L /* beyond ANSI C, see man feature_test_macros */
#include <netdb.h> /* getaddrinfo etc. */
#undef _POSIX_C_SOURCE

#include <sys/socket.h>
#include <sys/un.h> /* UNIX sockets */
#include <sys/mman.h> /* memory management */
#include <netinet/tcp.h> /* TCP_NODELAY */
#include <fcntl.h> /* O_RDWR */
#include <unistd.h> /* sleep */

#define INVALID_SOCKET -1
#define SOCKET_ERROR -1

static int fd;
#endif

#include "elphy.h"

void sread(void *data, const int len) {
    int all, new;

    for (all = 0; all < len; all += new)
        if ((new = recv(fd, (char *) data + all, len - all, 0)) == SOCKET_ERROR)
            error("Cannot read from socket.");
}

void swrite(const void *data, const int len) {
    if (send(fd, data, len, 0) == SOCKET_ERROR)
        error("Cannot write to socket.");
}

void open_inet_socket(const char *host, const char *port) {
    const int yes = 1;
    struct addrinfo hints = {0}, *res, *r;

    #ifdef _WIN32
    WORD version;
    WSADATA wsadata;

    version = MAKEWORD(1, 0);

    if (WSAStartup(version, &wsadata))
        error("Cannot initiate use of Winsock DLL");
    #endif

    hints.ai_family = AF_UNSPEC; /* IPv4 or IPv6 */
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    hints.ai_protocol = 0; /* any protocol */

    if (getaddrinfo(host, port, &hints, &res))
        error("Cannot get address info.");

    for (r = res;; r = r->ai_next ? r->ai_next : (sleep(1), res))
        if ((fd = socket(r->ai_family, r->ai_socktype, r->ai_protocol))
                != INVALID_SOCKET) {
            if (!setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, (void *) &yes,
                    sizeof yes))
                if (!connect(fd, r->ai_addr, r->ai_addrlen)) {
                    freeaddrinfo(res);
                    return;
                }
            close(fd);
        }
}

void close_socket(void) {
    close(fd);
    #ifdef _WIN32
    WSACleanup();
    #endif
}

#ifndef _WIN32
void open_unix_socket(const char *host, const char *prefix) {
    struct sockaddr_un addr = {0};

    addr.sun_family = AF_UNIX;
    strncat(addr.sun_path, prefix, sizeof addr.sun_path - 1);
    strncat(addr.sun_path, host, sizeof addr.sun_path - strlen(prefix) - 1);

    if ((fd = socket(AF_UNIX, SOCK_STREAM, 0)) == -1)
        error("Cannot create UNIX socket.");

    while (connect(fd, (struct sockaddr *) &addr, sizeof addr))
        sleep(1);
}

void *shm_attach(const int len) {
    char name[256] = "/";
    int namelen, mfd;
    void *addr;

    sread(&namelen, sizeof namelen);
    sread(name + 1, namelen);
    name[namelen + 1] = '\0';

    if ((mfd = shm_open(name, O_RDWR, 0)) == -1)
        error("Cannot open shared memory");

    if ((addr = mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_SHARED, mfd, 0))
            == MAP_FAILED)
        error("Cannot map shared memory");

    close(mfd);

    return addr;
}

void shm_detach(void *addr, const int len) {
    munmap(addr, len);
}
#endif
