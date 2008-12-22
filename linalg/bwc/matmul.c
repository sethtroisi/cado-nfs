#include "matmul.h"

/* these are just trampoline functions. A macro in the .h would serve the
 * same purpose, but could be misleading. Here at least, ctags bring the
 * reader quickly to the fact that we have macros.
 */

matmul_ptr matmul_build(abobj_ptr x, const char * filename)
{
    return MATMUL(build)(x, filename);
}

matmul_ptr matmul_reload_cache(abobj_ptr x, const char * filename)
{
    return MATMUL(reload_cache)(x, filename);
}

void matmul_save_cache(matmul_ptr mm, const char * filename)
{
    MATMUL(save_cache)(mm, filename);
}

void matmul_mul(matmul_ptr mm, abt * dst, abt const * src, int d)
{
    MATMUL(mul)(mm, dst, src, d);
}

void matmul_report(matmul_ptr mm) { MATMUL(report)(mm); }
void matmul_clear(matmul_ptr mm) { MATMUL(clear)(mm); }
