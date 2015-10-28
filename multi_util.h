#ifndef MULTI_UTIL_H
#define MULTI_UTIL_H

/*
 *   multi_util.h  version 12
 *    -- procedures shared among components
 */
#include "maf.h"

struct str_node {
    char* str;
    struct str_node* next;
};

void do_revcompl(char *s, int len);
void rev_comp(struct mafComp* c, int textSize);
void rc(struct mafAli *a);


// keep alignment starting from position beg on reference
struct mafAli* keep_ali(struct mafAli* ali, int beg);

// print part of an alignment ali from beg to end which are
// on the first component, according to col number
int print_part_ali_col(struct mafAli* ali, int beg, int end, FILE* fp);
struct mafAli* make_part_ali_col(struct mafAli* ali, int cbeg, int cend);

// locate column of struct mafAli with a given position in sequence 1
int mafPos2Col(struct mafComp *c, int pos, int textSize);

/* ------------ start of code to manage the output list of alignments ----------
*  The output alignments are not guarateed to be found in increasing order of start
*  position in the reference sequence.  As they are generated, we buffer them into a
*  properly sorted list, and frequently flush the list.
*/

struct mafAli* retrieve_first(struct mafAli** head);
void seperate_cp_wk(struct mafAli** cp_list, struct mafAli** wk_list, char* chr);

int colPos2Maf_after(struct mafComp* comp, int col);
int colPos2Maf_before(struct mafComp* comp, int col);

void parseSrcName(char* srcName, char* name, char* src);
void parseSrcName2(struct mafComp*);

int overlap(int beg1, int end1, int beg2, int end2);

#endif
