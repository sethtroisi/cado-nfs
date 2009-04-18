/* Reads a file with n 64-bit vectors (each one encoded as a
   16-digit hexadecimal string).
   Outputs a transposed file with 64 n-bit vectors (each one encoded as
   ceil(n/64) 16-digit hexadecimal strings) */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <inttypes.h>
#include "utils.h"

unsigned int read_hexstring(FILE *f, uint64_t * ptr, unsigned int n)
{
    unsigned int i = 0;
    unsigned int ir = 0;
    uint64_t x = (uint64_t) 0UL;
    uint64_t v;
    char c;

    for( ; isspace(c = fgetc(f)) ; );
    ungetc(c, f);

    for( ; i < n ; ) {
        char c;
        c = fgetc(f);

        if (c >= '0' && c <= '9') v = c - '0';
        else if (c >= 'a' && c <= 'f') v = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v = c - 'A' + 10;
        else {
            ungetc(c, f);
            return i;
        }

        /* 0 <= v < 16 */
        x |= v << ir;
        i += 4;
        ir += 4;
        if (ir == 64) {
            ir = 0;
            *ptr++ = x;
            x = (uint64_t) 0UL;
        }
    }
    if (i == n)
        return n;

    /* now n < i < n + 4. According to the specification, we check that
     * the top i-n bits are zero, and return n.
     */
    v >>= n+4-i;
    ASSERT_ALWAYS(!v);
    return n;
}

int main(int argc, char * argv[])
{
    FILE * f;
    int i, j;
    unsigned long k;
    uint64_t * x;
    unsigned int size = 0;
    unsigned int alloc = 0;
    uint64_t w;


    /* print command line */
    fprintf (stderr, "# %s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    /* read all data in memory */
    alloc = 1024; x = malloc(alloc * sizeof(uint64_t));
    if (argc != 2) {
        fprintf(stderr, "missing input file\n");
        exit(1);
    }
    f = fopen(argv[1], "r");
    ASSERT_ALWAYS(f);
    fprintf(stderr, "reading %s...", argv[1]); fflush(stderr);

    for ( ; !(feof(f) || ferror(f)) ; size++) {
        if (size >= alloc) {
            /* fprintf(stderr, "realloc(%d)", alloc); fflush(stderr); */
            alloc <<= 1;
            x = realloc(x, alloc * sizeof(uint64_t));
        }
        read_hexstring (f, x + size, 64);
        if (feof(f) || ferror(f))
            break;
    }

    fprintf(stderr, "read %d rows\n", size);

    fclose(f);

    for(i = 0 ; i < 64 ; i++) {
        j = 0;
        w = 0;
        for(k = 0 ; k < size ; k++) {
            w |= (uint64_t) ((x[k] >> i) & 1UL) << j;
            j ++;
            if (j == 64) {
                printf ("%016"PRIx64" ", w);
                j = 0;
                w = 0;
            }
        }
        /* print last word */
        if (j > 0)
            printf ("%" PRIx64, w);
        printf("\n");
    }

    fprintf (stderr, "Printed %d bits per dependency\n", size);
    return 0;
}
