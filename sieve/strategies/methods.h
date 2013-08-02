#ifndef __METHODS_H__
#define __METHODS_H__

#include <stdio.h>

#define METHOD_MAXBITS 50

//////////////////////////////////////////////////////////////////////////////
// The type for histograms of probabilities.
// Typically, the i-th coefficient is the probability related to primes
// with i-bits.

typedef float histogram_t [METHOD_MAXBITS+1];

//////////////////////////////////////////////////////////////////////////////
// Defines for the congruences modulo 12

#define MOD12_1  0
#define MOD12_5  1
#define MOD12_7  2
#define MOD12_11 3

//////////////////////////////////////////////////////////////////////////////
// Method types
// ECM, p-1, p+1
// FIXME: this should not be duplicated from facul.h !!!

#define PM1 0
#define PP1_27 1
#define PP1_65 2
#define ECM11 3
#define ECMM12 4

//////////////////////////////////////////////////////////////////////////////
// Main struct for a method
typedef struct {
    int type;
    int B1;
    int B2;
    int param;              // could be sigma
    float ms[3];            // millisecs, for i words input
    histogram_t success[4]; // proba to find a prime for each congruence mod 12
} cofac_method_struct;

typedef cofac_method_struct  cofac_method_t[1];
typedef cofac_method_struct *cofac_method_ptr;
typedef const cofac_method_struct *cofac_method_srcptr;


//////////////////////////////////////////////////////////////////////////////
// Protos

int method_read_stream(cofac_method_ptr meth, FILE *f);

// return number of methods read from file f, and stored in *meths (array
// allocated by the function, must be free'd by caller).
int methods_read(cofac_method_t **meths, const char *f);

void method_print(cofac_method_srcptr meth, FILE *file);
void method_print_full(cofac_method_srcptr meth, FILE *file);

#endif   /* __METHODS_H__ */
