#include "maf.h"
#include "util.h"
#include "multi_util.h"
#include "align_util.h"

static int MAX_SPECIES = 1000;
int CONNECTION_THRESHOLD = 50;
int SAME_CONNECTION = 30;

struct uAli** aliList2uAliArr(struct mafAli* subtreeAli, int treeCount) {
    struct mafAli *ali, *next;
    char* unused;
    struct uAli** subtree;
    int i, textSize, j;

    if ( treeCount == 0 )
        return NULL;

    subtree = (struct uAli**)malloc(treeCount*sizeof(struct uAli*));
    subtree[0] = (struct uAli*)malloc(treeCount*sizeof(struct uAli));
    for (i=0, ali=subtreeAli; i<treeCount; i++, ali=next) {
        next = ali->next;
        ali->next = NULL;
        subtree[i] = subtree[0] + i;
        subtree[i]->ali = ali;
        subtree[i]->index = i;
        subtree[i]->flipped = 'n';
        subtree[i]->sortContig = NULL;
        subtree[i]->next = NULL;
        subtree[i]->oali = NULL;
        subtree[i]->indexes = NULL;
        textSize = ali->textSize;
        unused = subtree[i]->usedArr = (char*)malloc(textSize*sizeof(char));
        for (j=0; j<textSize; j++)
            unused[j] = 'u';
    }
    return subtree;
}


int compar_uAli_start(const void* a, const void* b) {
    return (*((struct uAli**)a))->start - (*((struct uAli**)b))->start;
}

int compar_uAli_src(const void* a, const void* b) {
    return strcmp((*((struct uAli**)a))->sortContig, (*((struct uAli**)b))->sortContig);
}

int sort_uAli_contigs(struct uAli** uAliArr, int arrSize) {
    int i, front;
    char* prevStr;

    if ( arrSize == 0 )
        return 0;

    qsort(uAliArr, arrSize, sizeof(struct uAli*), compar_uAli_src);
    for (i=1, front=0, prevStr=uAliArr[0]->sortContig; i<arrSize; i++) {
        if ( strcmp(uAliArr[i]->sortContig, prevStr) != 0 ) { // end of a contig
            prevStr = uAliArr[i]->sortContig;
            qsort(uAliArr+front, i-front, sizeof(struct uAli*), compar_uAli_start);
            front = i;
        }
    }
    qsort(uAliArr+front, i-front, sizeof(struct uAli*), compar_uAli_start);  // the last segment
    return 0;
}


struct sortuAlis* do_sortuAlis(struct uAli** uAliArr, int totalSize, char* name, int subSize) {
    struct sortuAlis* suAlis;
    struct mafComp* comp;
    int i, k;

    suAlis = (struct sortuAlis*)malloc(sizeof(struct sortuAlis));
    if ( subSize > 0 )
        suAlis->sortuAliArr = (struct uAli**)malloc(subSize*sizeof(struct uAli*));
    else
        suAlis->sortuAliArr = NULL;
    suAlis->sortuAliArrSize = subSize;
    suAlis->sortedSpeciesName = copy_string(name);

    for (i=k=0; i<totalSize; i++) {
        for (comp = uAliArr[i]->ali->components; comp!=NULL; comp=comp->next)
            if ( strcmp(comp->name, name)==0 )
                break;
        if ( comp != NULL ) {
            suAlis->sortuAliArr[k++] = uAliArr[i];
            if ( uAliArr[i]->sortContig != NULL ) {
                free( uAliArr[i]->sortContig );
                uAliArr[i]->sortContig = NULL;
            }
            uAliArr[i]->sortContig = copy_string(comp->contig);
            if ( comp->strand == '+' )
                uAliArr[i]->start = comp->start;
            else
                uAliArr[i]->start = comp->srcSize - comp->start - comp->size;
            uAliArr[i]->end = uAliArr[i]->start + comp->size - 1;
        }
    }
    if ( k != subSize )
        fatalf("Inconsistent sizes k: %d and subSize: %d\n", k, subSize);
    sort_uAli_contigs(suAlis->sortuAliArr, subSize);

    suAlis->fronts = (int*)malloc(subSize*sizeof(int));
    suAlis->ends   = (int*)malloc(subSize*sizeof(int));
    for (k=0; k<subSize; k++) {
        suAlis->fronts[k] = suAlis->sortuAliArr[k]->start;
        suAlis->ends[k] = suAlis->sortuAliArr[k]->end;
    }
    return suAlis;
}

void initialize_uAliFile(struct mafAli* head, struct uAliFile* tAli) {
    char* names[MAX_SPECIES];
    struct mafComp* comp;
    struct mafAli* ali;
    int speciesCount=0, uAliCount=0, i;

    for (ali=head; ali != NULL; ali = ali->next) {
        for (comp=ali->components; comp != NULL; comp=comp->next) {
            for (i=0; i<speciesCount && i<1000; i++)
                if ( strcmp(comp->name, names[i]) == 0 )
                    break;
            if ( i == speciesCount )
                names[speciesCount++] = copy_string(comp->name);
        }
        uAliCount++;
    }

    tAli->speciesCount = speciesCount;
    if ( speciesCount > 0 ) {
        tAli->speciesAliCount = (int*)malloc(speciesCount*sizeof(int));
        tAli->speciesNames = (char**)malloc(speciesCount*sizeof(char*));
        tAli->sorted = (struct sortuAlis**)malloc(speciesCount*sizeof(struct sortuAlis*));
    } else {
        tAli->speciesAliCount = NULL;
        tAli->speciesNames = NULL;
        tAli->sorted = NULL;
    }
    tAli->uAliCount = uAliCount;

    for (i=0; i<speciesCount; i++) {
        (tAli->speciesNames)[i] = names[i];
        (tAli->speciesAliCount)[i] = 0;
    }

    for (ali=head; ali!=NULL; ali=ali->next)
        for (comp=ali->components; comp!=NULL; comp=comp->next) {
            for (i=0; i<speciesCount; i++)
                if ( strcmp(comp->name, names[i])==0 )
                    break;
            if ( i==speciesCount )
                fatalf("non-included species: %d\n", i);
            (tAli->speciesAliCount)[i]++;
        }

    tAli->uAliArr = aliList2uAliArr(head, uAliCount);

    for (i=0; i<speciesCount; i++)
        (tAli->sorted)[i] = do_sortuAlis(tAli->uAliArr, uAliCount, (tAli->speciesNames)[i], (tAli->speciesAliCount)[i]);
}

/* Given: intervals cbeg1-cend1 and cbegN-cendN, where
*    (1) positions cbeg1-cend1 are contained in the top row of alignment leftali, must be positive orient
*    (3) positions cbegN-cendN are contained in the top row of alignment rightali. might be negative orient
*/
int connectionAgreement2(struct mafAli* leftali, struct mafAli* rightali, int cbeg1, int cend1, int cbegN, int cendN, struct pwuAliFiles* pws) {
    struct mafAli* pw;
    struct mafComp *compA, *compB, *compa=NULL, *compb=NULL;
    char *topspecies, *botspecies;
    struct position a, b, c, d;
    int i, j, k, pairK=pws->pairK, marker1, marker2, overbeg, overend, cbeg, cend, beg1, end1, beg2, end2, leftK, rightK, tmp;
    int *existConnections, expectConnection, existConnection=0;
    int ab_mid_y, cd_mid_y, overmid;

    if ( leftali->components->strand == '-' )
        fatalf("left top component is not positive orientation: %s\n", leftali->components->name);

    for (leftK=0, compA=leftali->components; compA != NULL; compA=compA->next)
        leftK++;
    for (rightK=0, compB=rightali->components; compB != NULL; compB=compB->next)
        rightK++;

    existConnections = (int*)malloc(pairK*sizeof(int));
    for (i=0; i<pairK; i++)
        existConnections[i] = 0;
    expectConnection = leftK*rightK;

    for (compA=leftali->components; compA!=NULL; compA=compA->next) {
        marker1 = 0;
        if ( compA->strand == '-' ) {
            rev_comp(compA, leftali->textSize);
            for (compB=rightali->components; compB!=NULL; compB=compB->next)
                rev_comp(compB, rightali->textSize);
            tmp = cendN;
            cendN = rightali->textSize - cbegN - 1;
            cbegN = rightali->textSize - tmp - 1;
            tmp = cend1;
            cend1 = leftali->textSize - cbeg1 -1;
            cbeg1 = leftali->textSize - tmp - 1;
            marker1 = 1;
        }
        for (compB=rightali->components; compB!=NULL; compB=compB->next) {
            for (i=0; i<pairK; i++) {
                if ( pws->pwuAliFileArrs[i]->uAliCount == 0 )
                    continue;
                if ( pws->pwuAliFileArrs[i]->speciesCount < 2 )
                    fatal("pairwise alignment species number less than 2\n");
                topspecies = pws->pwuAliFileArrs[i]->speciesNames[0];
                botspecies = pws->pwuAliFileArrs[i]->speciesNames[1];
                if ( ( strcmp(compA->name, topspecies) == 0 && strcmp(compB->name, botspecies) == 0 )
                        || ( strcmp(compA->name, botspecies) == 0 && strcmp(compB->name, topspecies) == 0 ) )
                    break;
            }
            if ( i == pairK)
                continue;

            for (k=0; k<pws->pwuAliFileArrs[i]->speciesCount; k++)
                if ( strcmp(compA->name, pws->pwuAliFileArrs[i]->sorted[k]->sortedSpeciesName) == 0 )
                    break;
            if ( k == pws->pwuAliFileArrs[i]->speciesCount )
                fatalf("no sorted species: %s\n", compA->name);

            for (j=0; j<pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArrSize; j++) {
                if ( pws->pwuAliFileArrs[i]->sorted[k]->fronts[j] > (compA->start + compA->size - 1 ) )
                    continue;
                if ( pws->pwuAliFileArrs[i]->sorted[k]->ends[j] < compA->start )
                    continue;
                pw = pws->pwuAliFileArrs[i]->sorted[k]->sortuAliArr[j]->ali;
                if ( strcmp(compA->name, pw->components->name)==0 ) {
                    compa = pw->components;
                    compb = compa->next;
                } else {
                    compb = pw->components;
                    compa = compb->next;
                }
                if ( strcmp(compa->contig, compA->contig) != 0 || strcmp(compb->contig, compB->contig)!=0 )
                    continue;
                if ( compa->strand == '+')
                    if ( compb->strand != compB->strand )
                        continue;
                marker2 = 0;
                if ( compa->strand == '-') {
                    if ( compb->strand == compB->strand )
                        continue;
                    rev_comp(compa, pw->textSize);
                    rev_comp(compb, pw->textSize);
                    marker2 = 1;
                }

                a.x = beg2 = colPos2Maf_after(compA, cbeg1);
                b.x = end2 = colPos2Maf_before(compA, cend1);

                overbeg = (beg2 > compa->start ? beg2 : compa->start );
                overend = (end2 < compa->start+compa->size-1 ? end2 : compa->start+compa->size-1);
                if ( overbeg > overend )
                    continue;

                a.y = beg1 = colPos2Maf_after(compB, cbegN);
                b.y = end1 = colPos2Maf_before(compB, cendN);

                cbeg = mafPos2Col(compa, overbeg, pw->textSize);
                cend = mafPos2Col(compa, overend, pw->textSize);
                beg2 = colPos2Maf_after(compb, cbeg);
                end2 = colPos2Maf_before(compb, cend);

                if ( overlap(beg1, end1, beg2, end2) == 1 ) {
                    c.x = compa->start;
                    c.y = compb->start;
                    d.x = compa->start + compa->size - 1;
                    d.y = compb->start + compb->size - 1;
                    overbeg = ( a.x > c.x ? a.x : c.x );
                    overend = ( b.x < d.x ? b.x : d.x );
                    overmid = (overbeg + overend)/2;
                    ab_mid_y = b.y - (b.x - overmid)*(b.y - a.y)/(double)(b.x - a.x);
                    cd_mid_y = d.y - (d.x - overmid)*(d.y - c.y)/(double)(d.x - c.x);
                    if ( ab_mid_y - cd_mid_y >= -SAME_CONNECTION && ab_mid_y - cd_mid_y <= SAME_CONNECTION )
                        existConnections[i] = 1;
                }
                if ( marker2 == 1 ) {
                    rev_comp(compa, pw->textSize);
                    rev_comp(compb, pw->textSize);
                }
            }
        } // for compB
        if ( marker1 == 1 ) {
            rev_comp(compA, leftali->textSize);
            for (compB=rightali->components; compB!=NULL; compB=compB->next)
                rev_comp(compB, rightali->textSize);
            tmp = cendN;
            cendN = rightali->textSize - cbegN - 1;
            cbegN = rightali->textSize - tmp - 1;
            tmp = cend1;
            cend1 = leftali->textSize - cbeg1 -1;
            cbeg1 = leftali->textSize - tmp - 1;
        } // if (marker1==1)
    } // for compA

    for (i=0; i<pairK; i++)
        if ( existConnections[i] == 1 )
            existConnection++;

    if ( existConnection*100/expectConnection >= CONNECTION_THRESHOLD )
        return 1;
    return 0;
}
