#include <sys/types.h>
#include <stdio.h>
#include "macros.h"
#include "types.h"
#include "auxfuncs.h"
#include "variables.h"

int computed_m_param;
int computed_n_param;
int total_work;
coord_t nrows;


void consistency_check(const char * name, int computed, int hardcoded)
{
	if (computed==hardcoded)
		return;

	die("The variable %s is hard-coded "
		" to the value %d.\n"
		"The parameters you supplied imply a value of %d. The program"
		" fails\n",1,name,hardcoded,computed);
}
