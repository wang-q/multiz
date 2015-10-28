#ifndef MZ_SCORES_H
#define MZ_SCORES_H

#include "maf.h"

int	**ss, // 128x128 array of substitution scores
*gop; // 16-position array of quasi-natural gap-open penalties

int gap_open, gap_extend;

#define SS(c,d) ss[c][d]
#define GAP(s,t,u,v) gop[(s<<3)+(t<<2)+(u<<1)+v]
#define GAP2(s,t,u,v) GAP((s == '-'), (t == '-'), (u == '-'), (v == '-'))

void init_scores70();

double mafScoreRange(struct mafAli *maf, int start, int size);

#endif
