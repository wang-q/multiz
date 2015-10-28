/*
 *  multi_util.c version 12
*/

#include <string.h>
#include "util.h"
#include "maf.h"
#include "mz_scores.h"
#include "multi_util.h"

int radius = 30;
int MIN_OUTPUT_WID = 1;
int LRG_BREAK_WID = 20;
int SML_BREAK_WID = 2;

int OVERLAP_THRESHOLD = 50;
int MIN_CHAIN = 300;
int MIN_CLUSTER_CHAIN = 300;
int OVERLAP_LEN_THREH = 300;
int MIN_DISTANCE = 400;
int MIN_SPB = -100;
int row2 = 0;

int force = 0;
int execute = 1;
int verbose = 1;

char *PREFIX = NULL;
char *OPERAT = NULL;
char *USER_PATH = NULL;

// return a mafAli structure based on the input ali starting on position beg
struct mafAli *keep_ali(struct mafAli *ali, int beg) {
    int len, col_beg, count, i;
    char *s;
    struct mafComp *comp, *del_ptr;

    len = strlen(ali->components->text);

    col_beg = mafPos2Col(ali->components, beg, ali->textSize);
    for (; col_beg > 0 && ali->components->text[col_beg - 1] == '-'; col_beg--);

    for (comp = ali->components; comp != NULL;) {
        for (count = i = 0; i < col_beg; i++)
            if (comp->text[i] != '-')
                count++;                                 // number of bases before beg
        if (comp->size - count < 1) {             // delete
            if (comp != ali->components) {
                for (del_ptr = ali->components;
                     del_ptr->next != NULL && del_ptr->next != comp; del_ptr = del_ptr->next);
                if (del_ptr != NULL)
                    del_ptr->next = comp->next;
                mafCompFree(&comp);
                comp = del_ptr->next;
            } else {
                ali->components = comp->next;
                mafCompFree(&comp);
                comp = ali->components;
            }
            continue;
        }
        comp->start = comp->start + count;          // modified on Aug. 12th, cut starting beg, (count exclu des beg)
        comp->size = comp->size - count;
        s = (char *) malloc((len - col_beg + 2) * sizeof(char));
        for (i = col_beg; i < len; i++)
            s[i - col_beg] = comp->text[i];   // col_beg records pos before beg
        s[len - col_beg] = '\0';
        free(comp->text);
        comp->text = s;
        comp = comp->next;
    }
    ali->textSize = len - col_beg;
    ali->score = mafScoreRange(ali, 0, len - col_beg);
    return ali;
}


struct mafAli *make_part_ali_col(struct mafAli *ali, int cbeg, int cend) {
    struct mafComp *comp, *ncomp, *pcomp;
    int beg, cols, chs, i;
    struct mafAli *nali;

    if (cend - cbeg + 1 == 0)
        return NULL;
    nali = (struct mafAli *) malloc(sizeof(struct mafAli));
    nali->components = NULL;
    nali->next = NULL;
    nali->textSize = cend - cbeg + 1;
    nali->score = mafScoreRange(ali, cbeg, cend - cbeg + 1);

    for (comp = ali->components; comp != NULL; comp = comp->next) {
        for (beg = comp->start - 1, cols = 0; cols < cbeg; cols++)
            if (comp->text[cols] != '-')
                beg++;
        beg++;
        for (cols = cbeg, chs = 0; cols <= cend; cols++)
            if (comp->text[cols] != '-')
                chs++;
        if (chs == 0)
            continue;    // no all-dashs rows
        ncomp = mafCpyComp(comp);
        ncomp->start = beg;
        ncomp->size = chs;
        ncomp->text = (char *) malloc((cend - cbeg + 2) * sizeof(char));
        for (i = cbeg; i <= cend; i++)
            ncomp->text[i - cbeg] = comp->text[i];
        ncomp->text[cend - cbeg + 1] = '\0';
        if (nali->components == NULL)
            nali->components = ncomp;
        else {
            for (pcomp = nali->components; pcomp->next != NULL; pcomp = pcomp->next);
            pcomp->next = ncomp;
        }
    }
    if (nali->components != NULL) {
        nali = mafColDashRm(nali);
        if (nali != NULL)
            nali->score = mafScoreRange(nali, 0, nali->textSize);
        return nali;
    } else {
        mafAliFree(&nali);
        return NULL;
    }
}

// beg end are on cols, starting at 0
int print_part_ali_col(struct mafAli *ali, int cbeg, int cend, FILE *fp) {
    struct mafAli *nali;

    nali = make_part_ali_col(ali, cbeg, cend);
    if (nali != NULL) if (row2 == 0 || nali->components->next != NULL)
        mafWrite(fp, nali);
    mafAliFree(&nali);
    return 0;
}


// pos starts at 0, col starts at 0
int mafPos2Col(struct mafComp *c, int pos, int textSize) {
    int col, p;

    if (pos < c->start || pos >= c->start + c->size)
        fatalf("mafPos2Col: %d not in %d-%d",
               pos, c->start, c->start + c->size - 1);
    for (col = 0, p = c->start - 1; col < textSize; ++col)
        if (c->text[col] != '-' && ++p == pos)
            break;

    return col;
}

struct mafAli *retrieve_first(struct mafAli **head) {
    struct mafAli *ali;

    if (*head == NULL)
        return NULL;
    ali = *head;
    *head = (*head)->next;
    ali->next = NULL;
    return ali;
}

void seperate_cp_wk(struct mafAli **cp_list, struct mafAli **wk_list, char *chr) {
    struct mafAli *prev, *last = NULL, *a, *b;

    for (prev = a = *cp_list; a != NULL;) {
        if (strcmp(chr, a->components->src) == 0) {//move to wk_list
            if (a == *cp_list) { // at head
                *cp_list = (*cp_list)->next;
                a->next = NULL;
                b = a;
                prev = a = *cp_list;
            } else {
                prev->next = a->next;
                a->next = NULL;
                b = a;
                a = prev->next;
            }
            if (*wk_list == NULL)
                last = *wk_list = b;
            else {
                last->next = b;
                last = b;
            }
        } else {
            prev = a;
            a = a->next;
        }
    }
}

// name and src are allocated arrays, srcName is to be parsed
void parseSrcName(char *srcName, char *name, char *src) {
    char *ptr;
    int len;

    for (ptr = srcName; *ptr != '\0' && *ptr != '.'; ptr++);

    len = ptr - srcName;
    strncpy(name, srcName, len);
    name[len] = '\0';
    if (*ptr == '\0' || *(ptr + 1) == '\0')
        strcpy(src, name);
        //    fatal("srcName in wrong format 1");
    else {
        ++ptr;
        strcpy(src, ptr);
    }
}

// name and src are not allocated
void parseSrcName2(struct mafComp *c) {
    char *ptr, bk;

    for (ptr = c->src; *ptr != '\0' && *ptr != '.'; ptr++);
    bk = *ptr;
    *ptr = '\0';
    c->name = copy_string(c->src);
    *ptr = bk;
    if (*ptr == '\0' || *(ptr + 1) == '\0')
        c->contig = copy_string(c->src);
    else {
        ++ptr;
        c->contig = copy_string(ptr);
    }
}


