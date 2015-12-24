#ifndef MAF_H
#define MAF_H
/*
 *   maf.h version 12
 */
// mini-maf.h - Multiple alignment format, header file stolen from Jim Kent and abused
#include <stdio.h>
#include <zlib.h>

#define MAX_INT ((int)(~(1<<(sizeof(int)*8-1))))
#define MIN_INT ((int)(1<<(sizeof(int)*8-1)))

/* A file full of multiple alignments. */
struct mafFile {
    struct mafFile *next;
    int version;        // Required
    char *scoring;        // Name of scoring scheme.
    struct mafAli *alignments;    // Possibly empty list of alignments.
    char *fileName;
    int line_nbr;
    int verbose;
    gzFile *fp;            // Open file if any. NULL except while parsing.
};

void mafFileFree(struct mafFile **pObj);

/* Free up a maf file including closing file handle if necessary. */

/* A multiple alignment. */
struct mafAli {
    struct mafAli *next;
    double score;
    struct mafComp *components;
    /* List of components of alignment */
    int textSize;
    /* Size of text in each component. */
    int chain_len;
};

void mafAliFree(struct mafAli **pObj);

/* Free up a maf alignment. */

/* A component of a multiple alignment. */
struct mafComp {
    struct mafComp *next;
    char *name;
    /* comman name of sequence source. */
    char *src;
    /* Name of sequence source.  */
    char *text;
    /* The sequence including dashes. */
    char *contig;
    int *mafPosMap;
    int srcSize;
    /* Size of sequence source.  */
    int start;
    /* Start within sequence. Zero based. If strand is - is relative to src end. */
    int size;
    /* Size in sequence (does not include dashes).  */
    short nameID;
    char strand;
    /* Strand of sequence.  Either + or -*/
    char paralog;
};

void mafCompFree(struct mafComp **pObj);

/* Free up a maf component. */

struct mafFile *mafOpen(char *fileName, int verbose);

/* Open up a .maf file for reading.  Read header and
 * verify. Prepare for subsequent calls to mafNext().
 * Prints error message and aborts if there's a problem. */

struct mafAli *mafNext(struct mafFile *mafFile);

/* Return next alignment in file or NULL if at end.
 * This will close the open file handle at end as well. */

struct mafFile *mafReadAll(char *fileName, int verbose);

/* Read in full maf file */

void mafWriteStart(FILE *f, char *scoring);

/* Write maf header and scoring scheme name (may be null) */

void mafWriteEnd(FILE *f);

void mafWrite(FILE *f, struct mafAli *maf);

/* Write next alignment to file. */

struct mafAli *mafColDashRm(struct mafAli *a);

struct mafComp *mafCpyComp(struct mafComp *);

#endif
