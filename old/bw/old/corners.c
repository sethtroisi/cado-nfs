#include <stdlib.h>
#include <string.h>
#include "macros.h"
#include "options.h"
#include "corners.h"

int process_corner(char ** argv,
		struct extended_option_desc * conf,
		void *p_val)
{
	char *endptr;
	((corner_s *)p_val)->i	= strtol(*argv,&endptr,10);
	if (*endptr!=',' && *endptr!=':') return FALSE;
	endptr++;
	((corner_s *)p_val)->j	= strtol(endptr,&endptr,10);
	return strlen(endptr)==0;
}
