#ifndef UTIL_H
#define UTIL_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>     /* HUGE_VAL */
#include <limits.h>   /* INT_MAX, INT_MIN, LONG_MAX, LONG_MIN, etc. */
#include <assert.h>
#include <errno.h>

typedef unsigned char uchar;

extern char *argv0;

void print_argv0(void);

#ifdef __GNUC__     /* avoid some "foo might be used uninitialized" warnings */

void fatal(const char *msg) __attribute__ ((noreturn));

void fatalf(const char *fmt, ...) __attribute__ ((noreturn));

void fatalfr(const char *fmt, ...) __attribute__ ((noreturn));

#else
void fatal(const char *msg);
void fatalf(const char *fmt, ...);
void fatalfr(const char *fmt, ...);
#endif

FILE *ckopen(const char *name, const char *mode);

void *ckalloc(size_t amount);

char *copy_string(const char *s);

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))

#endif
