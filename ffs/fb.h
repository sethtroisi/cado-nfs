#ifndef __FB_H__
#define __FB_H__

#include "types.h"

typedef struct {
    unsigned int size;  // nb of entries in the factor base
    fbideal_t *elts;
} factorbase_t;

// Read a factor base from file. Return 1 iff success.
// If maxdeg!=0, when reaching an ideal of degree > maxdeg, stop reading
// and return.
int fbread(factorbase_t * FB, const char *file, int maxdeg);


#endif   /* __FB_H__ */
