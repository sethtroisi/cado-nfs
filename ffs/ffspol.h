#ifndef __FFSPOL_H__
#define __FFSPOL_H__

#include "fppol.h"

typedef struct {
    int deg;
    int alloc;
    fppol_t *coeffs;
} ffspol_t;


#endif   /* __FFSPOL_H__ */
