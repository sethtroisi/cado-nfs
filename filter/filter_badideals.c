#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "portability.h"
#include "filter_common.h"
#include "filter_badideals.h"
#include "mod_ul.h"


#define MAXLINE 1024

void read_bad_ideals_info(const char *filename, allbad_info_t info)
{
    FILE *file = fopen(filename, "r");
    ASSERT_ALWAYS(file != NULL);

    info->n = 0;
    info->badid_info = NULL;

    char str[1024], *ptr, *nptr;
    while (fgets(str, MAXLINE, file)) {
        if (feof(file))
            break;
        if (str[0] == '#' || str[0] == '\n')
            continue;
        errno = 0;
        p_r_values_t p = strtoul(str, &ptr, 10);
        ASSERT_ALWAYS(errno == 0);
        unsigned long k = strtoul(ptr, &nptr, 10);
        ASSERT_ALWAYS(errno == 0);
        p_r_values_t rk = strtoul(nptr, &ptr, 10);
        ASSERT_ALWAYS(errno == 0);

        badid_info_struct_t item;
        item.p = p;
        item.k = k;
        item.rk = rk;
        item.pk = p;
        for (unsigned int i = 1; i < k; ++i)
            item.pk *= p;
        if (rk < item.pk)
            item.r = item.rk % p;
        else {
            p_r_values_t x = item.rk-item.pk;
            item.r = p + (x % p);
        }
        item.ncol = 0;
        do {
            long v = strtol(ptr, &nptr, 10);
            ASSERT_ALWAYS(errno == 0);
            if (ptr == nptr)
                break;
            item.val[item.ncol] = v;
            item.ncol++;
            ptr = nptr;
        } while (1);

        info->n++;
        info->badid_info = (badid_info_struct_t *) realloc(
                info->badid_info,
                (info->n)*sizeof(badid_info_struct_t));
        info->badid_info[info->n-1] = item;
    }
    fclose(file);
}

static unsigned long divmod(long a, long b, unsigned long m)
{
    int sig=0;
    if (a<0) {
        sig = 1;
        a = -a;
    }
    if (b<0) {
        sig = 1-sig;
        b = -b;
    }
    modulusul_t mm;
    residueul_t aa, bb;
    modul_initmod_ul(mm, m);
    modul_init(aa, mm); modul_init(bb, mm);
    modul_set_ul(aa, a, mm);
    modul_set_ul(bb, b, mm);
    modul_inv(bb, bb, mm);
    modul_mul(aa, aa, bb, mm);
    a = modul_get_ul(aa, mm);
    modul_clear(aa, mm); modul_clear(bb, mm);
    modul_clearmod(mm);
    if (sig) {
        a = (m-a);
        if (a == (long)m)
            a = 0;
    }
    return a;
}


void
handle_bad_ideals (MAYBE_UNUSED int *exp_above, int64_t a, uint64_t b,
        unsigned long p, MAYBE_UNUSED int e, allbad_info_t info)
{
    p_r_values_t r;
    if ((b % p) == 0)
        r = p+divmod(b,a,p);
    else
        r = divmod(a,b,p);
    for(int i = 0; i < info->n; ++i) {
        if (p != info->badid_info[i].p)
            continue;
        if (r != info->badid_info[i].r)
            continue;
        p_r_values_t rk;
        p_r_values_t pk = info->badid_info[i].pk;
        if ((b % p) == 0)
            rk = pk+divmod(b,a,pk);
        else
            rk = divmod(a,b,pk);
        if (rk != info->badid_info[i].rk)
            continue;
        for (unsigned int j = 0; j < info->badid_info[i].ncol; ++j) {
            int v = info->badid_info[i].val[j];
            if (v>=0)
                exp_above[j] = v;
            else {
                ASSERT_ALWAYS(e >= -v);
                exp_above[j] = e+v;
            }
        }
        return;
    }
    ASSERT_ALWAYS(0);
}

