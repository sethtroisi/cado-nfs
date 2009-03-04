#define _GNU_SOURCE     // sometimes we use elctric_alloc->mmap->MMAP_ANONYMOUS

#include <stdio.h>
#include "matmul.h"

#include "abase.h"
#include "macros.h"

int main(int argc, char * argv[])
{
    abobj_t xx MAYBE_UNUSED;
    abobj_init(xx);

    matmul_t mm;

    if (argc != 2) {
        fprintf(stderr, "Usage: ./build <file>\n");
        exit(1);
    }

    mm = matmul_build(xx, argv[1]);
    matmul_save_cache(mm, argv[1]);
    matmul_clear(mm);
}

