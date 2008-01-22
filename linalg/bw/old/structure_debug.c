#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "lingen_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "field_def.h"
#include "field_usage.h"

#ifndef NDEBUG

extern int magma_display;
int magma_display=0;

#define INCLUDING_STRUCTURE_DEBUG_AUTOMATIC_H_
#include "structure_debug_automatic.h"

#endif	/* NDEBUG */
		
