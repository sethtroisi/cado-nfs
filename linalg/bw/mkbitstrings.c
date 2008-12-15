#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <inttypes.h>
#include "manu.h"
#include "utils.h"

unsigned int read_hexstring(FILE *f, unsigned long * ptr, unsigned int n)
{
    unsigned int i = 0;
    unsigned int ir = 0;
    unsigned long x = 0UL;
    unsigned long v;
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

        x |= v << ir;
        i += 4;
        ir += 4;
        if (ir == ULONG_BITS) {
            ir = 0;
            *ptr++ = x;
            x = 0UL;
        }
    }
    if (i == n)
        return n;

    /* now n < i < n + 4. According to the specification, we check that
     * the top i-n bits are zero, and return n.
     */
    v >>= n+4-i;
    BUG_ON(v);
    return n;
}

unsigned int write_hexstring(FILE * f, const unsigned long * ptr, unsigned int n)
{
    const char hxits[] = "0123456789abcdef";
    unsigned int i = 0;
    for( ; i < n ; ) {
        unsigned long v = *ptr++;
        unsigned int ir;
        for(ir = 0 ; ir < ULONG_BITS && i < n ; ) {
            fputc(hxits[v&0x0f], f);
            if (ferror(f) || feof(f)) return i;
            v >>= 4;
            ir += 4;
            i += 4;
        }
    }
    return n;
}

int main(int argc, char * argv[])
{
    FILE * f;
    int i, j;
    unsigned long k;
    unsigned long * x;
    unsigned int size = 0;
    unsigned int alloc = 0;
    uint64_t w;


    /* print command line */
    fprintf (stderr, "# %s.r%s", argv[0], VERSION);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    /* read all data in memory */
    alloc = 1024; x = malloc(alloc * sizeof(unsigned long));
    if (argc != 2) {
        fprintf(stderr, "missing input file\n");
        exit(1);
    }
    f = fopen(argv[1], "r");
    BUG_ON(f == NULL);
    fprintf(stderr, "reading %s...", argv[1]); fflush(stderr);

    for ( ; !(feof(f) || ferror(f)) ; size++) {
        if (size >= alloc) {
            /* fprintf(stderr, "realloc(%d)", alloc); fflush(stderr); */
            alloc <<= 1;
            x = realloc(x, alloc * sizeof(unsigned long));
        }
        read_hexstring(f, x + size, 64);
        if (feof(f) || ferror(f))
            break;
    }

    fprintf(stderr, "read %d rows\n", size);

    fclose(f);


    for(i = 0 ; i < 64 ; i++) {
        j = 0;
        w = 0;
        for(k = 0 ; k < size ; k++) {
            w |= ((x[k] >> i) & 1UL) << j;
            j ++;
            if (j == 64) {
                printf ("%016"PRIu64" ", w);
                j = 0;
                w = 0;
            }
        }
        /* print last word */
        if (j > 0)
            printf ("%" PRIu64, w);
        printf("\n");
    }

    fprintf (stderr, "Printed %d bits per dependency\n", size);
    return 0;
}
