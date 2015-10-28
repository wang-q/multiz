#ifndef MULTI_UTIL_H
#define MULTI_UTIL_H

/*
 *   multi_util.h  version 12
 *    -- procedures shared among components
 */
#include "maf.h"

// keep alignment starting from position beg on reference
struct mafAli *keep_ali(struct mafAli *ali, int beg);

// print part of an alignment ali from beg to end which are
// on the first component, according to col number
int print_part_ali_col(struct mafAli *ali, int beg, int end, FILE *fp);

struct mafAli *make_part_ali_col(struct mafAli *ali, int cbeg, int cend);

// locate column of struct mafAli with a given position in sequence 1
int mafPos2Col(struct mafComp *c, int pos, int textSize);

/* ------------ start of code to manage the output list of alignments ----------
*  The output alignments are not guarateed to be found in increasing order of start
*  position in the reference sequence.  As they are generated, we buffer them into a
*  properly sorted list, and frequently flush the list.
*/

struct mafAli *retrieve_first(struct mafAli **head);

void seperate_cp_wk(struct mafAli **cp_list, struct mafAli **wk_list, char *chr);

void parseSrcName(char *srcName, char *name, char *src);

void parseSrcName2(struct mafComp *);

#endif
