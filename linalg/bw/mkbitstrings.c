#include <stdio.h>
#include <string.h>
#include "manu.h"

unsigned int write_hexstring(FILE * f, const unsigned long * ptr, unsigned int n)
{
    const char hxits[] = "0123456789abcdef";
    unsigned int i = 0;
    for( ; i < n ; ) {
        unsigned long v = *ptr++;
        for(unsigned int ir = 0 ; ir < ULONG_BITS && i < n ; ) {
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
    for(int i = 0 ; i < argc ; i++) {
        ix[i] = fopen(argv[i], "r");
        BUG_ON(ix == NULL);
    }
    unsigned long foo[BITS_TO_WORDS(argc, ULONG_BITS)];
    for(;;) {
        memset(foo,0,sizeof(foo));
        int fail = 0;
        for(int i = 0 ; i < argc ; i++) {
            unsigned long x;
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
    for(int i = 0 ; i < argc ; i++) {
        fclose(ix[i]);
    }
    return 0;
}
