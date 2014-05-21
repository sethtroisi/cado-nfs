#include "cado.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <string.h>

#include "portability.h"
#include "macros.h"
#ifdef  HAVE_CURL
#include "balancing_curl_source.h"
#endif
#include "balancing_file_source.h"

/* This is a stupid program for checking a matrix-times-vector product:
 *
 * ./a.out <matrix file> <vector file> : compute M * v
 * ./a.out -t <matrix file> <vector file> : compute v * M
 *
 * dimensions are guessed from the vector size. data files are expected
 * in binary format, and printed in the same binary format to stdout.
 * Pipe through xxd -ps -c8 to get binary hex (or xxd -ps -c8 | uniq -c
 * tto verify that you reach 0 as expected for the final kernel vector).
 */

void short_matmul(FILE * out, FILE * v,  const char * uri, int mul_left)
{
    struct stat sbuf[1];
    if (fstat(fileno(v), sbuf) < 0) {
        perror("input vector");
        exit(1);
    }

    unsigned long ncols = sbuf->st_size;
    unsigned long nrows;
    assert(ncols % sizeof(uint64_t) == 0);
    ncols /= sizeof(uint64_t);
    if (mul_left == 0) {
        fprintf(stderr, "Found %lu columns in input vector\n", ncols);
        nrows = 0;
    } else {
        fprintf(stderr, "Found %lu rows in input vector\n", ncols);
        nrows = ncols;
        ncols += ncols / 10;
    }
    uint64_t * vec;
    vec = malloc(ncols * sizeof(uint64_t));
    memset(vec, 0, ncols * sizeof(uint64_t));
    assert(vec);

    data_source_ptr m;
#ifdef  HAVE_CURL
    int is_curl;
    if ((is_curl = (strstr(uri, "://") != NULL)))
        m = curl_source_alloc(uri, 0);
    else
#endif  /* HAVE_CURL */
    m = file_source_alloc(uri, 0);

    if (mul_left == 0) {
        unsigned long l = fread(vec, sizeof(uint64_t), ncols, v);
        if (l < ncols) {
            fprintf(stderr, "input vector: short read\n");
            exit(1);
        }
        fclose(v);
        for( ; ; ) {
            uint32_t x;
            uint32_t * ptr = &x;
            int k;
            k = m->get(m, &ptr, 1);
            ASSERT_ALWAYS(ptr == &x);
            uint32_t rowlen = x;
            if (k != 1) {
                if (k) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                break;
            }
            uint64_t w = 0;
            for(uint32_t i = 0 ; i < rowlen ; i++) {
                k = m->get(m, &ptr, 1);
                ASSERT_ALWAYS(ptr == &x);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                assert(x < ncols);
                w ^= vec[x];
            }
            // printf("%016" PRIx64 "\n", w);
            fwrite(&w, sizeof(uint64_t), 1, out);
        }
    } else {
        unsigned long rrows;
        for(rrows = 0 ; ; rrows++) {
            uint32_t x;
            uint32_t * ptr = &x;
            int k;
            k = m->get(m, &ptr, 1);
            uint32_t rowlen = x;
            if (k != 1) {
                if (k) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                break;
            }
            uint64_t w = 0;
            k = fread(&w, sizeof(uint64_t), 1, v);
            if (k != 1) {
                fprintf(stderr, "matrix: short read\n");
                exit(1);
            }
            for(uint32_t i = 0 ; i < rowlen ; i++) {
                k = m->get(m, &ptr, 1);
                ASSERT_ALWAYS(ptr == &x);
                if (k != 1) {
                    fprintf(stderr, "matrix: short read\n");
                    exit(1);
                }
                if (x >= ncols) {
                    unsigned long old_ncols = ncols;
                    for( ; x >= ncols ; ncols += ncols / 10);
                    fprintf(stderr, "re-allocating vec array, now %lu cols\n", ncols);
                    vec = realloc(vec, ncols * sizeof(uint64_t));
                    memset(vec + old_ncols, 0, (ncols - old_ncols) * sizeof(uint64_t));
                }
                vec[x] ^= w;
            }
        }
        assert(nrows == rrows);
        unsigned long k;
        for(k = 0 ; k < nrows ; k++) {
            uint64_t w = vec[k];
            // printf("%016" PRIx64 "\n", w);
            fwrite(&w, sizeof(uint64_t), 1, out);
        }
    }

#ifdef  HAVE_CURL
    if (is_curl)
        curl_source_free(m);
    else
#endif
    file_source_free(m);
    free(vec);
}

void usage()
{
    fprintf(stderr, "Usage: short-matmul <matrix file> <vector file>\n");
    exit(1);
}
int main(int argc, char * argv[])
{
    int mul_left = 0;
    for( ; argc > 3 ; ) {
        if (strcmp(argv[1], "-t") == 0) {
            mul_left = 1;
            argv++, argc--;
        }
    }
    if (argc != 3) {
        usage();
    }
    /*
    FILE * m = fopen(argv[1], "rb");
    if (m == NULL) {
        perror(argv[1]);
        exit(1);
    }
    */
    FILE * v = fopen(argv[2], "rb");
    if (v == NULL) {
        perror(argv[1]);
        exit(1);
    }
    
    short_matmul(stdout, v, argv[1], mul_left);
}

