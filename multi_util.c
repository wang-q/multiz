/*
 *  multi_util.c version 12
*/

#include <string.h>
#include "util.h"
#include "seq.h"
#include "maf.h"
#include "mz_scores.h"
#include "multi_util.h"
//#include "myDebug.h"

static const char rcsid[] = "$Id: multi_util.c 142 2008-11-12 18:55:23Z rico $";

int radius=30;
int MIN_OUTPUT_WID=1;
int LRG_BREAK_WID=20;
int SML_BREAK_WID=2;

int OVERLAP_THRESHOLD=50;
int MIN_CHAIN=300;
int MIN_CLUSTER_CHAIN=300;
int OVERLAP_LEN_THREH=300;
int MIN_DISTANCE=400;
int MIN_SPB=-100;
int row2=0;

int force=0;
int execute=1;
int verbose=1;

char* PREFIX=NULL;
char* OPERAT=NULL;
char* USER_PATH=NULL;

static char dna_compl[256] =
    "                                             -                  "
    " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
    "                                                                "
    "                                                                ";
/* .............................................-.................. */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */

void do_revcompl(char *s, int len) {
    char c, *p = s + len - 1;

    while (s<=p) {
        c = dna_compl[(uchar)*s];
        *s = dna_compl[(uchar)*p];
        *p = c;
        ++s, --p;
    }
}

void rev_comp(struct mafComp* c, int textSize) {
    c->start = c->srcSize - (c->start + c->size);
    c->strand = (c->strand == '-' ? '+' : '-');
    do_revcompl(c->text, textSize);
}

// replace every component by its reverse complement
void rc(struct mafAli *a) {
    struct mafComp *c;

    for (c = a->components; c != NULL; c = c->next)
        rev_comp(c, a->textSize);
}

/*
void rc_uAli(struct uAli *a) {
  int len, i, j;
  char used;

  rc(a->ali);
  len = a->ali->textSize;
  for (i=0, j=len-1; i<=j; i++, j--) {
    used = a->used[i];
    a->used[i] = a->used[j];
    a->used[j] = used;
  }
  a->flipped = ( a->flipped == 'n' ? 'y' : 'n' );
}
*/

int parse_fasta(char *line, char **ident, int *start, int *end, int *srcSize) {
    char *s, *t, *copy = copy_string(line);
    int ret;

    if (copy[0] != '>' || (s = strchr(copy, ':')) == NULL ||
            (t = strchr(s+1, ':')) == NULL ||
            sscanf(t+1, "%d-%d:%*c:%d", start, end, srcSize) != 3)
        ret = 0;
    else {
        *t = '\0';
        if (strchr(copy, ' ') != NULL)
            ret = 0;
        else {
            ret = 1;
            *ident = copy_string(s+1);
        }
    }
    free(copy);
    return ret;
}

char *seq_ident(char *filename) {
    FILE *fp = fopen(filename, "r");
    char buf[500], *ident, *retur;
    int x, y, z;

    if (fp == NULL || fgets(buf, 500, fp) == NULL ||
            parse_fasta(buf, &ident, &x, &y, &z) == 0)
        retur = filename;
    else
        retur = ident;
    if (fp != NULL)
        fclose(fp);
    return retur;
}

static int compar(const void *a, const void *b) {
    return (*((struct mafAli **)a))->components->start -
           (*((struct mafAli **)b))->components->start;
}

struct mafAli **mafArray(struct mafAli *a, int *nali) {
    struct mafAli *b, **A;
    int n, i;

    for (n = 0, b = a; b != NULL; b = b->next)
        ++n;
    *nali = n;
    A = ckalloc((n+10000)*sizeof(struct mafAli *));
    for (b = a, i = 0; i < n; ++i, b = b->next) {
        if (b == NULL)
            fatal("died in mafArray");
        A[i] = b;
    }
    if (b != NULL)
        fatal("did not reach end in mafArray");
    /* sort alignments by start of first component */
    qsort((void *)A, n, sizeof(struct mafAli *), compar);

    for (i = 1; i < n; ++i)
        if (A[i-1]->components->start > A[i]->components->start)
            fatal("failure in mafArray");

    return A;
}

void order_rows_species(struct mafAli **a, struct mafComp **location, char **species,
                        int nspecies) {
    struct mafComp *c, *last;
    char src[500], name[500];
    int i;

    for (i = 0; i < nspecies; ++i)
        location[i] = NULL;
    for (c = (*a)->components; c != NULL; c = c->next) {
        parseSrcName(c->src, name, src);
        for (i = 0; i < nspecies; ++i) {
            if (same_string(name, species[i]))
                break;
        }
        if (i == nspecies)
            continue;
        if (location[i] != NULL) {
            mafWrite(stderr, *a);
            fatalf("species `%s'named more than once",
                   c->src);
        }
        location[i] = c;
    }

    for (last = NULL, i = 0; i < nspecies; ++i)
        if (location[i] != NULL) {
            if (last == NULL)
                (*a)->components = location[i];
            else
                last->next = location[i];
            last = location[i];
        }
    if (last == NULL) {
        mafAliFree(a);
        *a = NULL;
    } else
        last->next = NULL;
}

/* not updated code -- 10/12/05
void order_rows(struct mafAli **a, struct mafComp **location, char **species,
		int nspecies, struct str_node* chrs, struct str_node* spes) {
	struct mafComp *c, *last;
	int i;
	struct str_node *chr, *spe;

	for (i = 0; i < nspecies; ++i)
		location[i] = NULL;
	for (c = (*a)->components; c != NULL; c = c->next) {
		for (i = 0; i < nspecies; ++i)
			if (same_string(c->name, species[i]) || same_string(c->src, species[i])) // !!!
				break;

		if (i == nspecies && chrs != NULL) {
		  for (chr=chrs, spe=spes; chr!=NULL; chr=chr->next, spe=spe->next)
		    if ( same_string(c->name, chr->str) || same_string(c->src, chr->str)) // !!!
		      break;
		  if ( spe != NULL) {
		    for (i=0; i<nspecies; i++)
		      if ( same_string(spe->str, species[i])) // !!!
			break;
		  }
		}
		if ( i==nspecies)
		  continue;

		if (location[i] != NULL) {
			mafWrite(stderr, *a);
			fatalf("species `%s'named more than once",
			       c->src);
		}
		location[i] = c;
	}

	for (last = NULL, i = 0; i < nspecies; ++i)
		if (location[i] != NULL) {
			if (last == NULL)
			  (*a)->components = location[i];
			else
				last->next = location[i];
			last = location[i];
		}
	if (last == NULL) {
	  mafAliFree(a);
	  *a = NULL;
	}
	else
	  last->next = NULL;
}
*/

/*
0 - ${COMMON_NAME}
1 - ${ENCODE_REGION}
2 - ${FREEZE_DATE}
3 - ${NCBI_TAXON_ID}
4 - ${ASSEMBLY_PROVIDER}
5 - ${ASSEMBLY_DATE}
6 - ${ASSEMBLY_ID}
7 - ${CHROMOSOME}
8 - ${CHROMOSOME_START}
9 - ${CHROMOSOME_END}
10 - ${CHROM_LENGTH}
11 - ${STRAND}
12 - ${ACCESSION}.${VERSION}
13 - ${NUM_BASES}
14 - ${NUM_N}
15 - ${THIS_CONTIG_NUM}
16 - ${TOTAL_NUM_CONTIGS}
*/
static int parseMSAHeader(char* fn, SEQ* sf, char* name, char* chr, int* start, int* end, char* strand, int* size) {
    char *arr[17], *ptr1, *ptr2;
    int i, len;

    for (ptr1=sf->header; *ptr1 != '>'; ptr1++)
        ;
    for (++ptr1, ptr2=ptr1, i=0; *ptr2 != '\0'; ptr2++) {
        if ( *ptr2 == '|') {
            len = ptr2 - ptr1;
            arr[i] = (char*)malloc((len+1)*sizeof(char));
            strncpy(arr[i], ptr1, len);
            arr[i][len] = '\0';
            i++;
            ptr1 = ptr2+1;
        }
    }
    if ( i == 17 ) {
        if ( strcmp(arr[0], ".")!=0 && strcmp(arr[7], ".")!=0 && strcmp(arr[8], ".")!=0 && strcmp(arr[9], ".")!=0 && strcmp(arr[10], ".")!=0 && strcmp(arr[11], ".")!=0 ) {
            strcpy(name, arr[0]);
            strcpy(chr,  arr[7]);
            *start = atoi(arr[8]);
            //*end   = atoi(arr[9]);
            *end = *start + strlen((char*)sf->seq) - 1;
            *size  = atoi(arr[10]);
            *strand = arr[11][0];
            //if ( strcmp(arr[6], ".") != 0 )
            //strcpy(name, arr[6]);
            return 0;
        }
        if ( strcmp(arr[0], ".")!=0 && strcmp(arr[15], ".")!=0 && strcmp(arr[13], ".")!=0 && strcmp(arr[11], ".")!=0) {
            strcpy(name, arr[0]);
            strcpy(chr,  arr[15]);
            *start = 1;
            *end   = atoi(arr[13]);
            *size  = *end;
            *strand = arr[11][0];
            return 0;
        }
    }
    strcpy(name, fn);
    strcpy(chr, fn);
    *start = 1;
    *size = *end = strlen((char*)sf->seq);
    *strand = '+';
    return 0;

    //fatalf("Wrong format header: %s", sf->header);
    //return 1;
}

int parseHeader(char* fn, SEQ* sf, char* name, char* chr, int* start, int* end, char* strand, int* size) {

    if ( sf->header[0] != '>' )
        fatal("Wrong fasta header!");
    if ( sscanf(sf->header, ">%[^:]:%[^:]:%d-%d:%c:%d", name, chr, start, end, strand, size) == 6 )
        return 0;
    if ( sscanf(sf->header, ">%[^:]:%[^:]:%d:%c:%d", name, chr, start, strand, size) == 5 ) {
        *end = *start + strlen((char*)sf->seq) - 1;
        return 0;
    } else
        return parseMSAHeader(fn, sf, name, chr, start, end, strand, size);
}

/*
// assume src is allocated enough space
int parseHeader(char* fn, SEQ* sf, char* name, char* chr, int* start, int* end, char* strand, int* size) {
  char buf[1000], *ptr, *p;
  int i, len, j, control;

  if ( sf->header[0] != '>' )
    fatal("Wrong fasta header.");

  i=0;
  control = 0;
  for (p=ptr=sf->header+1; ; ptr++) {
    if ( *ptr==':' || *ptr=='\0' ) {
      len = ptr-p;
      strncpy(buf, p, len);
      buf[len] = '\0';

      switch (i) {
      case 0: strcpy(name, buf);
	break;
      case 1: strcpy(chr, buf);
	break;
      case 2:
	for (j=0; j<len; j++)
	  if ( buf[j] < '0' || buf[j] > '9') {
	    control = 1;
	    break;
	  }
	if ( control == 1)
	  break;
	*start = atoi(buf);
	break;
      case 3:
	*strand = buf[0];
	if ( *strand != '+' && *strand != '-' )
	  control = 1;
	break;
      case 4:
	for (j=0; j<len; j++)
	  if ( buf[j] < '0' || buf[j] > '9') {
	    control = 1;
	    break;
	  }
	if ( control == 1)
	  break;
	*size = atoi(buf);
	break;
      default:
	break;
      }
      i++;
      p=ptr+1;
      if ( control == 1)
	break;
    }
    if ( *ptr == '\0')
      break;
  }
  if ( control == 0 && i==5) {
    *end = *start + strlen((char*)sf->seq) - 1;
    return 0;
  }

  return parseMSAHeader(fn, sf, name, chr, start, end, strand, size);
}
*/

/*
int splitContigs(char* seqfile, struct str_node** chrs, int* total) {
  char filename[500], *filestr, tmpBuf[500], name[200], chr[200], strand;
  int i, iter, start, end, size;
  SEQ* sf;
  FILE* fpw;
  struct str_node* n;

  *chrs = NULL;
  iter=*total=0;
  sf = seq_open(seqfile);
  if (!seq_read(sf))
    fatalf("Cannot read sequence from %s", seqfile);

  //assume sf->header is in format of ">name:postions", need to get between > :
  //splitting sequence file
  do {
    parseHeader(sf->header, name, chr, &start, &end, &strand, &size);
    sprintf(tmpBuf, "_Bz_%s", chr);


    for (i=0,n=*chrs; n!=NULL; n=n->next,i++) {
      if ( strcmp(tmpBuf, n->str)==0 ) {
	sprintf(filename, "%s.%d", tmpBuf, iter++);
	break;
      }
    }

    filestr = (char*)malloc(500*sizeof(char));
    if ( i<*total)
      strcpy(filestr, filename);
    else
      strcpy(filestr, tmpBuf);
    n = (struct str_node*)malloc(sizeof(struct str_node));
    n->str = copy_string(filestr);
    n->next = *chrs;
    *chrs = n;
    //    chrs[*total] = filestr;

    *total = *total + 1;
    fpw = ckopen(filestr, "w");
    fprintf(fpw, "%s\n", sf->header);
    fprintf(fpw, "%s", sf->seq);
    fclose(fpw);
  } while (seq_read(sf));
  seq_close(sf);
  return 0;
}
*/

void flip_comp(struct mafAli *a) {
    struct mafComp *c1, *c2;
    if ((c1 = a->components) == NULL)
        fatal("alignment does not have rows");
    if ((c2 = c1->next) == NULL)
        return;
    c1->next = c2->next;
    c2->next = c1;
    a->components = c2;
    if (c2->strand == '-') {
        c2->start = c2->srcSize - (c2->start + c2->size);
        c2->strand = '+';
        do_revcompl(c2->text, a->textSize);
        c1->start = c1->srcSize - (c1->start + c1->size);
        c1->strand = (c1->strand == '-' ? '+' : '-');
        do_revcompl(c1->text, a->textSize);
    }
}

// swap the first and second components of a struct mafAli
void flip_comps(struct mafAli *a) {
    for ( ; a != NULL; a = a->next) {
        flip_comp(a);
    }
}

// return a mafAli structure based on the input ali starting on position beg
struct mafAli* keep_ali(struct mafAli* ali, int beg) {
    int len, col_beg, count, i;
    char* s;
    struct mafComp* comp, *del_ptr;

    len = strlen(ali->components->text);

    col_beg = mafPos2Col(ali->components, beg, ali->textSize);
    for (; col_beg > 0 && ali->components->text[col_beg-1]=='-'; col_beg--)
        ;

    for (comp=ali->components; comp != NULL; ) {
        for (count=i=0; i<col_beg; i++)
            if ( comp->text[i] != '-')
                count++;                                 // number of bases before beg
        if ( comp->size - count < 1 ) {             // delete
            if ( comp != ali->components) {
                for ( del_ptr = ali->components; del_ptr->next != NULL && del_ptr->next != comp; del_ptr = del_ptr->next)
                    ;
                if ( del_ptr != NULL)
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
        s = (char*)malloc( (len-col_beg+2)*sizeof(char));
        for (i=col_beg; i<len; i++)
            s[i-col_beg] = comp->text[i];   // col_beg records pos before beg
        s[len-col_beg] = '\0';
        free(comp->text);
        comp->text = s;
        comp = comp->next;
    }
    ali->textSize = len-col_beg;
    ali->score = mafScoreRange(ali, 0, len-col_beg);
    return ali;
}


// beg, end are on the first component, not col number, starts at 0
int print_part_ali(struct mafAli* ali, int beg, int end, FILE* fp) {
    struct mafComp* comp, *ncomp, *last;
    int rel_beg, rel_end, cols, chs, col_beg, col_end, len, i;
    struct mafAli* nali;

    rel_beg = beg - ali->components->start;
    rel_end = end - ali->components->start;

    len = strlen(ali->components->text); // get col number for rel_beg
    for (cols=0,chs=col_beg=col_end=-1; cols < len && chs <= rel_end; cols++) {
        if ( ali->components->text[cols] != '-') {
            chs++;
            if ( chs==rel_beg)
                col_beg = cols;
            if ( chs==rel_end)
                col_end = cols;
        }
    }
    if ( col_beg < 0 || col_end < 0)
        return 0;

    nali = mafNewAli( mafScoreRange(ali, col_beg, col_end-col_beg+1), col_end-col_beg+1);
    last = NULL;

    for (comp=ali->components; comp != NULL; comp = comp->next) {
        for (beg=comp->start-1, cols=0; cols < col_beg; cols++)
            if ( comp->text[cols] != '-')
                beg++;
        beg++;
        for (cols = col_beg, chs=0; cols<=col_end; cols++)
            if ( comp->text[cols] != '-')
                chs++;
        if (chs==0) // all '-'
            continue;
        ncomp = mafCpyComp(comp);
        ncomp->start = beg;
        ncomp->size = chs;
        ncomp->text = (char*)malloc((col_end-col_beg+2)*sizeof(char));

        for (i=col_beg; i<=col_end; i++)
            ncomp->text[i-col_beg] = comp->text[i];
        ncomp->text[col_end-col_beg+1] = '\0';
        if ( nali->components == NULL)
            nali->components = last = ncomp;
        else {
            last->next = ncomp;
            last = ncomp;
        }
    }
    if ( nali->components != NULL)
        mafWrite(fp, nali);
    mafAliFree(&nali);
    return 0;
}

struct mafAli* make_part_ali_col(struct mafAli* ali, int cbeg, int cend) {
    struct mafComp* comp, *ncomp, *pcomp;
    int beg, cols, chs, i;
    struct mafAli* nali;

    if ( cend-cbeg+1 == 0 )
        return NULL;
    nali = (struct mafAli*)malloc(sizeof(struct mafAli));
    nali->components = NULL;
    nali->next = NULL;
    nali->textSize = cend-cbeg+1;
    nali->score = mafScoreRange(ali, cbeg, cend-cbeg+1);

    for (comp = ali->components; comp != NULL; comp=comp->next) {
        for (beg=comp->start-1, cols=0; cols < cbeg; cols++)
            if ( comp->text[cols] != '-')
                beg++;
        beg++;
        for (cols=cbeg, chs=0; cols<=cend; cols++)
            if ( comp->text[cols] != '-')
                chs++;
        if (chs==0)
            continue;    // no all-dashs rows
        ncomp = mafCpyComp(comp);
        ncomp->start = beg;
        ncomp->size = chs;
        ncomp->text = (char*)malloc((cend-cbeg+2)*sizeof(char));
        for (i=cbeg; i<=cend; i++)
            ncomp->text[i-cbeg] = comp->text[i];
        ncomp->text[cend-cbeg+1] = '\0';
        if ( nali->components==NULL)
            nali->components = ncomp;
        else {
            for (pcomp=nali->components; pcomp->next != NULL; pcomp = pcomp->next)
                ;
            pcomp->next = ncomp;
        }
    }
    if ( nali->components != NULL) {
        nali = mafColDashRm(nali);
        if ( nali != NULL)
            nali->score = mafScoreRange(nali, 0, nali->textSize);
        return nali;
    } else {
        mafAliFree(&nali);
        return NULL;
    }
}

// beg end are on cols, starting at 0
int print_part_ali_col(struct mafAli* ali, int cbeg, int cend, FILE* fp) {
    struct mafAli* nali;

    nali = make_part_ali_col(ali, cbeg, cend);
    if ( nali != NULL )
        if (row2==0 || nali->components->next != NULL)
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
    //if (p != pos)
    //fatal("mafPos2Col: cannot happen");
    return col;
}

int mafPos2Col_v2(struct mafComp *c, int pos) {
    int col, p;

    if (pos < c->start || pos >= c->start + c->size)
        fatalf("mafPos2Col_v2: %d not in %d-%d",
               pos, c->start, c->start + c->size - 1);
    for (col = 0, p = c->start - 1; c->text[col]!='\0'; ++col)
        if (c->text[col] != '-' && ++p == pos)
            break;
    //if (p != pos)
    //fatal("mafPos2Col: cannot happen");
    return col;
}

// Create an alignment equal to a segment of a given alignment; don't set score
struct mafAli *mafSlice(struct mafAli *a, int start_col, int beyond_col) {
    int i, len = beyond_col - start_col;
    struct mafAli *new = ckalloc(sizeof(struct mafAli));
    struct mafComp *ac, *nc, *prev = NULL;

    new->textSize = len;
    new->next = NULL;
    new->components = NULL;
    for (ac = a->components; ac != NULL; ac = ac->next) {
        nc = mafCpyComp(ac);
        for (i = 0; i < start_col; ++i)
            if (ac->text[i] != '-')
                ++(nc->start);
        nc->text = ckalloc((len+1)*sizeof(char));
        for (nc->size = i = 0; i < len; ++i)
            if ((nc->text[i] = ac->text[i+start_col]) != '-')
                ++(nc->size);
        nc->text[len] = '\0';
        if ( prev == NULL )
            new->components = nc;
        else
            prev->next = nc;
        prev = nc;
    }
    return new;
}

// flush output-list entries that start before a specified position
//struct mafAli *mafFlush(struct mafAli *output_list, int pos, FILE* fp) {
struct mafAli *mafFlush(struct mafAli *output_list, int pos) {
    struct mafAli *atmp;

    while (output_list != NULL && output_list->components->start <= pos) {
        mafWrite(stdout, output_list);
        atmp = output_list;
        output_list = output_list->next;
        mafAliFree(&atmp);
    }
    return output_list;
}

// insert an alignment into the output list, sorting by components->start
struct mafAli *mafInsert(struct mafAli *output_list, struct mafAli *a) {
    struct mafAli *m, *next;
    int beg1;

    if (a == NULL)
        return output_list;

    beg1 = a->components->start;
//fprintf(stderr, "mafInsert %d\n", beg1);
    if (output_list == NULL || output_list->components->start >= beg1) {
        // put a at the start of the output list
        a->next = output_list;
        return a;
    }
    for (m = output_list; (next = m->next) != NULL &&
            next->components->start < beg1;  m = next)
        ;
    // splice a between m and m->next
    a->next = next;
    m->next = a;
    return output_list;
}

void insert_ali(struct mafAli** ali_head, struct mafAli* ali) {
    ali->next = *ali_head;
    *ali_head = ali;
}

int compar_ali_top_start(const void* a, const void* b) {
    return (*(struct mafAli**)a)->components->start - (*(struct mafAli**)b)->components->start;
}

int change_neg_pos(struct mafComp* comp, int pos) { // pos is on comp
    if ( comp->strand == '+')
        return pos;
    return comp->srcSize-1-pos;
}

// currently break_ali cut ali with two rows
// bp is on top component maf position, ali2 is already allocated
void break_ali(struct mafAli* ali, int bp, struct mafAli* ali2) { // the second new ali includes bp if possible
    struct mafComp* comp, *srcComp;
    int col, new_start_col, new_end_col, new_start_top_maf, new_start_bot_maf, new_end_top_maf, new_end_bot_maf, top_maf, bot_maf;

    new_start_top_maf= new_start_bot_maf = new_end_top_maf = new_end_bot_maf =-1;
    col = mafPos2Col(ali->components, bp, ali->textSize);
    for (new_start_col = col; new_start_col<ali->textSize; new_start_col++ )
        if ( ali->components->text[new_start_col] != '-' && ali->components->next->text[new_start_col] != '-')
            break;
    for (new_end_col = col-1; new_end_col>=0; new_end_col--)
        if ( ali->components->text[new_end_col] != '-' && ali->components->next->text[new_end_col] != '-')
            break;
    top_maf=ali->components->start-1;
    bot_maf=ali->components->next->start-1;
    for (col=0; col<ali->textSize; col++) {
        if ( ali->components->text[col]!='-')
            top_maf++;
        if ( ali->components->next->text[col]!='-')
            bot_maf++;
        if ( col==new_end_col ) {
            new_end_top_maf = top_maf;
            new_end_bot_maf = bot_maf;
        }
        if ( col==new_start_col) {
            new_start_top_maf = top_maf;
            new_start_bot_maf = bot_maf;
            break;
        }
    }

    ali2->next = NULL;
    ali2->textSize = ali->textSize - new_start_col;

    srcComp = ali->components;
    comp = mafCpyComp(srcComp);
    comp->start = new_start_top_maf;
    comp->size = srcComp->size - (new_start_top_maf - srcComp->start);
    comp->text = (char*)malloc( (ali->textSize - new_start_col + 1)*sizeof(char));
    srcComp = srcComp->next;
    comp->next = mafCpyComp(srcComp);
    comp->start = new_start_bot_maf;
    comp->size = srcComp->size - (new_start_bot_maf - srcComp->start);
    comp->text = (char*)malloc( (ali->textSize - new_start_col+1)*sizeof(char));

    for (col=0; col<ali->textSize-new_start_col; col++) {
        comp->text[col] = ali->components->text[col+new_start_col];
        comp->next->text[col] = ali->components->next->text[col+new_start_col];
    }
    comp->text[ali->textSize-new_start_col] = comp->next->text[ali->textSize-new_start_col] = '\0';
    ali2->components = comp;
    ali2->score = mafScoreRange(ali2, 0, ali2->textSize);

    ali->components->size = new_end_top_maf - ali->components->start + 1;
    ali->components->next->size = new_end_bot_maf - ali->components->next->start + 1;
    ali->components->text[new_end_col+1] = '\0';
    ali->components->next->text[new_end_col+1] = '\0';

    ali->textSize = new_end_col+1;
    ali->score = mafScoreRange(ali, 0, ali->textSize);
}

struct mafAli* retrieve_first(struct mafAli** head) {
    struct mafAli* ali;

    if (*head == NULL)
        return NULL;
    ali = *head;
    *head = (*head)->next;
    ali->next = NULL;
    return ali;
}

void seperate_cp_wk(struct mafAli** cp_list, struct mafAli** wk_list, char* chr) {
    struct mafAli *prev, *last=NULL, *a, *b;

    for (prev=a=*cp_list; a!=NULL; ) {
        if ( strcmp(chr, a->components->src)==0) {//move to wk_list
            if (a==*cp_list) { // at head
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
            if ( *wk_list == NULL )
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

// if colth is '-', mafpos should be the first non - after that
int colPos2Maf_after(struct mafComp* comp, int col) {
    int i, pos;
    for ( pos = comp->start-1, i=0; i<col; i++)
        if (comp->text[i] != '-')
            pos++;
    pos++;
    if ( pos > comp->start + comp->size - 1)
        return -1;
    return pos;
}

// if colth is '-', mafpos should be the last non - before that
int colPos2Maf_before(struct mafComp* comp, int col) {
    int i, pos;
    for (pos = comp->start-1, i=0; i<=col; i++)
        if ( comp->text[i] != '-')
            pos++;
    if ( pos < comp->start)
        return -1;
    return pos;
}

int test_ali(struct mafAli* ali) {
    int i;
    char* ch;
    for ( ch=ali->components->text,i=0; *ch != '\0'; ch++)
        if ( *ch != '-' )
            i++;
    if ( i!=ali->components->size) {
        printf("ali->components->size: %d real-size: %d\n", ali->components->size, i);
        return 1;
    }
    for ( ch=ali->components->next->text, i=0; *ch!='\0'; ch++)
        if ( *ch != '-')
            i++;
    if ( i!=ali->components->next->size) {
        printf("ali->components->next->size: %d realsize: %d\n", ali->components->next->size, i);
        return 1;
    }
    return 0;
}

// name and src are allocated arrays, srcName is to be parsed
void parseSrcName(char* srcName, char* name, char* src) {
    char* ptr;
    int len;

    for (ptr=srcName; *ptr!='\0' && *ptr!='.'; ptr++)
        ;

    len = ptr - srcName;
    strncpy(name, srcName, len);
    name[len] = '\0';
    if ( *ptr == '\0' || *(ptr+1)=='\0' )
        strcpy(src, name);
    //    fatal("srcName in wrong format 1");
    else {
        ++ptr;
        strcpy(src, ptr);
    }
}

// name and src are not allocated
void parseSrcName2(struct mafComp* c ) {
    char* ptr, bk;

    for (ptr=c->src; *ptr!='\0' && *ptr!='.'; ptr++)
        ;
    bk = *ptr;
    *ptr = '\0';
    c->name = copy_string(c->src);
    *ptr = bk;
    if ( *ptr == '\0' || *(ptr+1)=='\0' )
        c->contig = copy_string(c->src);
    //    fatal("srcName in wrong format 1");
    else {
        ++ptr;
        c->contig = copy_string(ptr);
    }
}

int overlap(int beg1, int end1, int beg2, int end2) {
    int over_beg, over_end, over_len;
    double over_threshold;

    if ( beg2 > end1 || beg1 > end2 )
        return 0;

    over_beg = ( beg1 > beg2 ? beg1 : beg2);
    over_end = ( end1 < end2 ? end1 : end2);
    over_len = over_end - over_beg + 1;
    over_threshold = (double)OVERLAP_THRESHOLD/100;

    if ( (double)over_len/(end1-beg1+1) > over_threshold
            || (double)over_len/(end2-beg2+1) > over_threshold
            || over_len >= OVERLAP_LEN_THREH )
        return 1;
    return 0;
}


