#include "abase-u64.h"

#include "pad.h"
#define P(X)    PAD(abase_u64,X)
#define ABASE_F(t,n,a) t P(n) a

#include "abase-binary-dotprod-generic.c"
