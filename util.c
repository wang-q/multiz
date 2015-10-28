#include <stdarg.h>
#include "util.h"

char *argv0;

/* print_argv0 ---------------------------------------- print name of program */
void print_argv0(void) {
    if (argv0) {

        char *p = strrchr(argv0, '/');
        (void) fprintf(stderr, "%s: ", p ? p + 1 : argv0);
    }
}

/* fatal ---------------------------------------------- print message and die */
void fatal(const char *msg) {
    fatalf("%s", msg);
}

/* fatalf --------------------------------- format message, print it, and die */
void fatalf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    print_argv0();
    (void) vfprintf(stderr, fmt, ap);
    (void) fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

/* ckopen -------------------------------------- open file; check for success */
FILE *ckopen(const char *name, const char *mode) {
    FILE *fp;

    if ((fp = fopen(name, mode)) == NULL)
        fatalf("Cannot open %s.", name);
    return fp;
}

void *ckalloc(size_t amount) {
    void *p;

    if ((long) amount < 0) {                                  /* was "<= 0" -CR */
        fatal("ckalloc: request for negative space.");
    }
    if (amount == 0)
        amount = 1; /* ANSI portability hack */
    if ((p = malloc(amount)) == NULL)
        fatalf("Ran out of memory trying to allocate %lu.",
               (unsigned long) amount);
    return p;
}


/* copy_string ---------------------- save string s somewhere; return address */
char *copy_string(const char *s) {
    char *p = ckalloc(strlen(s) + 1);    /* +1 to hold '\0' */
    return strcpy(p, s);
    //strcpy(p, s);
    //return p;
}


void fatalfr(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);

    fflush(stdout);
    print_argv0();
    (void) vfprintf(stderr, fmt, ap);
    (void) fprintf(stderr, ": %s\n", strerror(errno));
    va_end(ap);
    exit(1);
}

