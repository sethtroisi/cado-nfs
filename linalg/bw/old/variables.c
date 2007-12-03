#include <sys/types.h>
#include <stdio.h>
#include "macros.h"
#include "types.h"
#include "auxfuncs.h"
#include "variables.h"

coord_t nrows;
coord_t ncols;
int computed_m_param;
int computed_n_param;
int bank_num;
int total_work;			/* # of iterations. */
int periodicity;			/* freq. of saves, roughly. */

void consistency_check(const char * name, int computed, int hardcoded)
{
	if (computed==hardcoded)
		return;

	die("The variable %s is hard-coded in h/slave/slave_params.h"
		" to the value %d.\n"
		"The parameters you supplied imply a value of %d. The program"
		" fails\n",1,name,hardcoded,computed);
}

