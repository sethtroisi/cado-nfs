#include "cado.h"

/*
#include "cado.h"
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/statvfs.h>


#include "bwc_config.h"
*/
#include "matmul.h"

/* Include all matmul implementations here */

#define MM_IMPL_MAGIC_basic     1
#define MM_IMPL_MAGIC_sliced    2
#define MM_IMPL_MAGIC_bucket    3
#define MM_IMPL_MAGIC_threaded  4

#define CPP_PAD2(A,B) A ## B
#define MM_IMPL_MY_MAGIC(X) CPP_PAD2(MM_IMPL_MAGIC_,X)

#if !defined(MM_IMPL)
#error "Please compile this file with the MM_IMPL macro defined"
#elif     MM_IMPL_MY_MAGIC(MM_IMPL) == MM_IMPL_MAGIC_basic
#include "matmul-basic.h"
#elif     MM_IMPL_MY_MAGIC(MM_IMPL) == MM_IMPL_MAGIC_sliced
#include "matmul-sliced.h"
#elif     MM_IMPL_MY_MAGIC(MM_IMPL) == MM_IMPL_MAGIC_bucket
#include "matmul-bucket.h"
#elif     MM_IMPL_MY_MAGIC(MM_IMPL) == MM_IMPL_MAGIC_threaded
#include "matmul-threaded.h"
#else
#error "Undefined MM implementation provided in MM_IMPL"
#endif

/* Now some cpp glue which sets up the different options */
#define MATMUL_NAME(kind,func) matmul_ ## kind ## _ ## func
#define MATMUL_T(func) matmul_ ## func ## _t

#define REBIND_F(mm, kind,func) \
        mm->bind->func = (MATMUL_T(func)) & MATMUL_NAME(kind, func)

#define SET_IMPL(mm, impl_) do { mm->bind->impl = # impl_ ; } while (0)

#define REBIND_ALL(mm, kind) do {					\
        REBIND_F(mm, kind, build_cache);				\
        REBIND_F(mm, kind, reload_cache);				\
        REBIND_F(mm, kind, save_cache);			        	\
        REBIND_F(mm, kind, mul);					\
        REBIND_F(mm, kind, report);					\
        REBIND_F(mm, kind, clear);					\
        REBIND_F(mm, kind, init);					\
        REBIND_F(mm, kind, auxv);					\
        REBIND_F(mm, kind, aux);					\
        SET_IMPL(mm, kind);                                             \
    } while (0)

void matmul_solib_do_rebinding(matmul_ptr mm)
{
    REBIND_ALL(mm, MM_IMPL);
}
