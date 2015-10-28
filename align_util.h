#ifndef ALIGN_UTIL_H
#define ALIGN_UTIL_H

/*
 *  align_util.h version 1
 *    -- structures and precedures shared by aligners
 */

#include "maf.h"

struct position {
    int x, y;
};

// transform a single ali to ordered uAli with a given list of species
struct ordereduAli {
    char** contigs;
    int* begins, *ends;  // positively oriented
};

//----< wrapper of ali structure with mark array >------
struct uAli {
    struct mafAli* ali; // ali contains textSize field
    char* usedArr;      // 'u' unused; 'o' occupied
    struct uAli* next;
    char* sortContig;
    struct ordereduAli *oali;
    int index, start, end, cbeg, cend, rows, *indexes, incoming;
    char flipped;
};

//----< sorted array of uAlis >-----
struct sortuAlis {
    struct uAli** sortuAliArr;
    int *fronts, *ends;
    char* sortedSpeciesName;
    int sortuAliArrSize;
};

//----< collection of uAlis from an alignment file >-----
struct uAliFile {
    struct uAli** uAliArr;
    struct sortuAlis** sorted;
    char** speciesNames;
    int*   speciesAliCount;
    char* filename;
    int uAliCount, speciesCount;
};

//----< collection of pariwise alignment files >------
struct pwuAliFiles {
    struct uAliFile** pwuAliFileArrs;
    int pairK;  // number of pairwise alignment files
};


struct uAli** aliList2uAliArr(struct mafAli* subtreeAli, int treeCount);

void initialize_uAliFile(struct mafAli* head, struct uAliFile* tAli);

int connectionAgreement2(struct mafAli* a2, struct mafAli* a3, int cbeg2, int cend2, int cbeg3, int cend3, struct pwuAliFiles* pws);

int sort_uAli_contigs(struct uAli** uAliArr, int arrSize);

struct sortuAlis* do_sortuAlis(struct uAli** uAliArr, int totalSize, char* name, int subSize);

#endif
