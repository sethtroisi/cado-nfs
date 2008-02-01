#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "manu.h"

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
    argc--, argv++;
    FILE * ix[argc];
    int i, j, nbits = 0;
    unsigned long x;

    for(i = 0 ; i < argc ; i++) {
        ix[i] = fopen(argv[i], "r");
        BUG_ON(ix == NULL);
    }
    uint64_t w;

#if 0
    unsigned long foo[BITS_TO_WORDS(argc, ULONG_BITS)];
    for(;;) {
        memset(foo,0,sizeof(foo));
        int fail = 0;
        for(i = 0 ; i < argc ; i++) {
            fscanf(ix[i], "%lu", &x);
            fail += feof(ix[i]) || ferror(ix[i]);
            foo[i / ULONG_BITS] ^= x << (i % ULONG_BITS);
        }
        if (fail == argc)
            break;
        BUG_ON(fail);
        write_hexstring(stdout, foo, argc);
        fputc('\n', stdout);
    }
#else
    /* format is one line per dependency */
    for (i = 0; i < argc; i++)
      {
        j = 0;
        w = 0;
        nbits = 0;
        while (!(feof(ix[i]) || ferror(ix[i])))
          {
            fscanf(ix[i], "%lu", &x);
            if (feof(ix[i]) || ferror(ix[i]))
              break;
            w |= x << j;
            j ++;
            nbits ++;
            if (j == 64)
              {
                printf ("%lx ", w);
                j = 0;
                w = 0;
              }
          }
        /* print last word */
        if (j > 0)
          printf ("%lx\n", w);
      }
#endif
    for(i = 0 ; i < argc ; i++) {
        fclose(ix[i]);
    }
    fprintf (stderr, "Printed %d bits per dependency\n", nbits);
    return 0;
}
