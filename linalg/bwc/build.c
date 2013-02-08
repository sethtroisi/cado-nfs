#include "cado.h"

#include <string.h>
#include <stdio.h>

#include "bwc_config.h"
#include "matmul.h"
#include "portability.h"
#include "macros.h"
#include "mpfq/mpfq.h"
#include "mpfq/abase_vbase.h"
#include "matmul-mf.h"

void usage()
{
    fprintf(stderr, "Usage: ./build <file> [<impl>] [right|left]\n");
    exit(1);
}
int main(int argc, char * argv[])
{
    abase_vbase xx;
    abase_vbase_oo_field_init_byfeatures(xx,
                /* MPFQ_PRIME_MPZ, (2 as mpz), */ /* this is by default I
                                                 think */
                MPFQ_GROUPSIZE, 64,
                MPFQ_DONE);

    matmul_t mm;

    const char * impl = "basic";

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    if (argc != 2 && argc != 3 && argc != 4) {
        usage();
    }

    if (argc >= 3)
        impl = argv[2];

    int d = 0;
    int nd = 2;

    if (argc == 4) {
        if (strcmp(argv[3],"left") == 0)  {
            d = 0;
        } else if(strcmp(argv[3],"right") == 0)  {
            d = 1;
        } else {
            usage();
        }
        nd = 1;
    }

    for( ; nd-- ; d ^= 1) {
        if (d == 1) {
            fprintf(stderr, "Saving cache for matrix-times-vector\n");
        } else {
            fprintf(stderr, "Saving cache for vector-times-matrix\n");
        }
        mm = matmul_init(xx, 0, 0, argv[1], impl, NULL, d);

        ASSERT_ALWAYS(mm->store_transposed == !d);
        matrix_u32 m;
        /* This excludes the with-coefficients case for now */
        mf_prepare_matrix_u32(mm, m, argv[1], 0);
        matmul_build_cache(mm, m);
        matmul_save_cache(mm);
        matmul_clear(mm);
    }
}
