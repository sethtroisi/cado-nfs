#include <string.h>
#include "las-report-stats.h"

void las_report_init(las_report_ptr p)
{
    memset(p, 0, sizeof(las_report));
}

void las_report_clear(las_report_ptr p)
{
    memset(p, 0, sizeof(las_report));
}

void las_report_accumulate(las_report_ptr p, las_report_ptr q)
{
    p->reports += q->reports;
    p->survivors0 += q->survivors0;
    p->survivors1 += q->survivors1;
    p->survivors2 += q->survivors2;
    p->ttsm    += q->ttsm;
    p->ttf     += q->ttf;
    for(int side = 0 ; side < 2 ; side++) {
        p->tn[side]  += q->tn[side];
    }
    for (int i1 = 0; i1 < 256; i1++) {
        for (int i2 = 0; i2 < 256; i2++) {
            p->survivor_sizes[i1][i2] += q->survivor_sizes[i1][i2];
            p->report_sizes[i1][i2] += q->report_sizes[i1][i2];
        }
    }
    memset(q, 0, sizeof(las_report));
}
