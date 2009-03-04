#include "abase-common.h"
#include "abase-u64n.h"

#include "pad.h"
#define P(X)    PAD(abase_u64n,X)
#define ABASE_F(t,n,a) t P(n) a

#include "abase-binary-dotprod-generic.c"
