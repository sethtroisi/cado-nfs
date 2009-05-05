#define _GNU_SOURCE     // sometimes we use elctric_alloc->mmap->MMAP_ANONYMOUS

#include <stdio.h>
#include "bwc_config.h"
#include "matmul.h"
#include "abase.h"
#include "macros.h"

int main(int argc, char * argv[])
{
    abobj_t xx MAYBE_UNUSED;
    abobj_init(xx);

    matmul_t mm;

    const char * impl = "basic";

    if (argc != 2 && argc != 3) {
        fprintf(stderr, "Usage: ./build <file> [<impl>]\n");
        exit(1);
    }

    if (argc == 3)
        impl = argv[2];

    fprintf(stderr, "Saving cache for matrix-times-vector\n");
    mm = matmul_build(xx, argv[1], argv[2], NULL, 1);
    fprintf(stderr, "Saving cache for vector-times-matrix\n");
    mm = matmul_build(xx, argv[1], argv[2], NULL, 0);

    matmul_save_cache(mm, argv[1]);
    matmul_clear(mm);
}

