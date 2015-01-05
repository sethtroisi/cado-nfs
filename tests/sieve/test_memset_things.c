#include "cado.h"
#include "las-norms.h"

int main()
{
    extern size_t max_cache, min_stos;
    max_cache = direct_write_vs_stos ();
    min_stos = stos_vs_write128 ();
    printf("movaps / rep-stosq cutoff: %zu(0x%zx) ;", min_stos, min_stos);
    printf(" rep-stosq / movntps cutoff:");
    if (max_cache != ~(size_t)0)
        printf(" %zu(0x%zx)\n", max_cache, max_cache);
    else
        printf(" never\n");
    return 0;
}

