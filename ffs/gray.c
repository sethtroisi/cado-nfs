#include<stdio.h>
#include<stdlib.h>
#include <assert.h>

/*
 Playing with Gray code. 
 Yes, gcc converts the __builtin_ctz() into the bsf instruction of the 
 x86 abi that is fast on recent Intel but is a complete sucker on AMD.
 (see Agner Fog's manual).
*/


// Small Gray code that can be unrolled if __builtin_ctz() is costly.
void printgk(unsigned int k)
{
    for (unsigned int i = 0; i < (1U<<k) - 1; ++i) {
        printf("%d\n", __builtin_ctz(i+1));
    }
}


// Print the successive xors, assuming we have already a function for 
// a Gray code of order k.
void print_incremental_gray(unsigned int n, unsigned int k)
{
    for (unsigned int i = 0; i < 1U<<(n-k); ++i) {
        printgk(k);
        printf("%d\n", k+__builtin_ctz(i+1));
    }
}


void binary_print(char * S, unsigned int x, int fix)
{
    unsigned int mask = 1U<<(fix-1);
    while (mask != 0) {
        if (x & mask) 
            S[0] = '1';
        else
            S[0] = '0';
        S++;
        mask >>= 1;
    }
    *S='\0';
}

unsigned int graycode(unsigned int n)
{
    return n ^ (n>>1);
}

#ifdef GRAY_MAIN

int main(int argc, char **argv)
{
    char S[1024], T[1024];
    if (argc != 2) {
        fprintf(stderr, "usage: ./a.out n k\n");
        return 1;
    }
    int K = atoi(argv[1]);
    int k = atoi(argv[2]);
    for (unsigned int i = 0; i < 1U<<K; ++i) {
        binary_print(S, i, K);
        binary_print(T, graycode(i), K);
        int offset = __builtin_ctz(i+1);
        assert (graycode(i+1) == (graycode(i) ^ (1U << offset)));
        printf("%u %s %s %d\n", i, S, T, offset);
    }
    printf("\n");
    print_incremental_gray(K, k);
    return 0;
}

#endif
