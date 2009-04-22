#include <stdarg.h>

#include "matmul.h"

#define MATMUL_DEFAULT_IMPL sliced

/* Include all matmul implementations here */
#include "matmul-basic.h"
#include "matmul-sliced.h"
#include "matmul-threaded.h"


/* Now some cpp glue which sets up the different options */
#define MATMUL_NAME(kind,func) matmul_ ## kind ## _ ## func
#define MATMUL_T(func) matmul_ ## func ## _t

#define REBIND_F(mm, kind,func) \
        mm->bind->func = (MATMUL_T(func)) & MATMUL_NAME(kind, func)
#define REBIND_ALL(mm, kind) do {					\
        REBIND_F(mm, kind, build);					\
        REBIND_F(mm, kind, reload_cache);				\
        REBIND_F(mm, kind, save_cache);			        	\
        REBIND_F(mm, kind, mul);					\
        REBIND_F(mm, kind, report);					\
        REBIND_F(mm, kind, clear);					\
        REBIND_F(mm, kind, auxv);					\
        REBIND_F(mm, kind, aux);					\
    } while (0)

/* these are just trampoline functions. A macro in the .h would serve the
 * same purpose, but could be misleading. Here at least, ctags bring the
 * reader quickly to the fact that we have macros.
 */

void do_rebinding(matmul_ptr mm, const char * impl)
{
#define CHECK_REBIND(mm, K) \
    if (strcmp(impl, #K) == 0) { REBIND_ALL(mm, K); } else

    if (impl == NULL) { REBIND_ALL(mm, MATMUL_DEFAULT_IMPL); } else
    CHECK_REBIND(mm, sliced)        // no semicolon !
    CHECK_REBIND(mm, threaded)      // no semicolon !
    CHECK_REBIND(mm, basic)         // no semicolon !
    { fprintf(stderr, "no implementation %s known\n", impl); exit(1); }
}

matmul_ptr matmul_build(abobj_ptr x, const char * filename, const char * impl, param_list pl)
{
    struct matmul_public_s fake[1];
    do_rebinding(fake, impl);
    matmul_ptr mm = fake->bind->build(x, filename, pl);
    if (mm == NULL) return NULL;
    do_rebinding(mm, impl);
    return mm;
}

matmul_ptr matmul_reload_cache(abobj_ptr x, const char * filename, const char * impl, param_list pl)
{
    struct matmul_public_s fake[1];
    do_rebinding(fake, impl);
    matmul_ptr mm = fake->bind->reload_cache(x, filename, pl);
    if (mm == NULL) return NULL;
    do_rebinding(mm, impl);
    return mm;
}

void matmul_save_cache(matmul_ptr mm, const char * filename)
{
    mm->bind->save_cache(mm, filename);
}

void matmul_mul(matmul_ptr mm, abt * dst, abt const * src, int d)
{
    mm->bind->mul(mm, dst, src, d);
}

void matmul_report(matmul_ptr mm) { mm->bind->report(mm); }
void matmul_clear(matmul_ptr mm) { mm->bind->clear(mm); }

void matmul_auxv(matmul_ptr mm, int op, va_list ap)
{
    mm->bind->auxv (mm, op, ap);
}

void matmul_aux(matmul_ptr mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_auxv (mm, op, ap);
    va_end(ap);
}
