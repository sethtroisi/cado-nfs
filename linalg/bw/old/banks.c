#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "params.h"
#include <assert.h>
#include "macros.h"
#include "types.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "filenames.h"
#include "variables.h"

char modulusname[MODULUSNAME_MAX_LENGTH];

static int w_get_header(char string[], char * filename)
{
	FILE	*f;
	int	res=0;

	f=fopen(filename,"r");
	if (f==NULL)
		return -1;
	if (fgets(string,DATA_IDENT_LEN,f)==NULL)
		res=-1;
	if (fclose(f)<0)
		return -1;
	else
		return res;
}

int try_to_write_bank(FILE ** p_fm, FILE ** p_fv)
{

#define fnm i_matrix_filename
#define fnv i_vector_filename
#define fnr i_result_filename

	if (exist(fnm)) {
		printf("Skipping bank %d (%s already there)\n",
				current_bank,fnm);
		return -1;
	} else if (exist(fnv)) {
		printf("Skipping bank %d (%s already there)\n",
				current_bank,fnv);
		return -1;
	} else if (exist(fnr)) {
		printf("Skipping bank %d (%s already there)\n",
				current_bank,fnr);
		return -1;
	}

	*p_fm = fopen(fnm,"w");
	*p_fv = fopen(fnv,"w");

#undef fnm
#undef fnv
#undef fnr

	return 0;
}
	
int try_to_load_bank(FILE ** p_fm, FILE ** p_fv)
{
	char idm[DATA_IDENT_LEN];
	char idv[DATA_IDENT_LEN];
	char *s,*t;
	int ex_m,ex_v,ex_r;
	int stub;
	int idlm,idlv;

	if (exist(wip_tag_filename)) {
		printf("Skipping bank %d (Work already in progress)\n"
		       "You may want to remove %s\n",current_bank,wip_tag_filename);
		return -1;
	}
	
#define fnm i_matrix_filename
#define fnv i_vector_filename
#define fnr i_result_filename

	ex_m=exist(fnm);
	ex_v=exist(fnv);
	ex_r=exist(fnr);

	if (! (ex_m && ex_v)) {
		printf("Skipping bank %d (incomplete data)\n",current_bank);
		return -1;
	}

	if (ex_r) {
		printf("Skipping bank %d (%s already there)\n",current_bank,fnr);
		return -1;
	}

	if (w_get_header(idm,fnm)<0 || w_get_header(idv,fnv)<0) {
		printf("Skipping bank %d (%s)\n",current_bank,strerror(errno));
		return -1;
	}

	idlm=strlen(idm);
	idlv=strlen(idv);


	printf("Bank %d : Identification strings :\n",current_bank);
	fputs(idm,stdout);
	fputs(idv,stdout);

	if (idlm!=idlv || strcmp(idm,idv)!=0) {
		printf("Skipping bank %d (IDs unmatched)\n",current_bank);
		return -1;
	}

	/* Get the matrix size */

	stub	= strcspn(idm,"0123456789");
	s	= idm + stub;

	nrows	= strtol(s,&t,10);
	stub	= strcspn(t,"0123456789");
	s	= t + stub;
	ncols	= strtol(s,&t,10);

	if (nrows == 0 || ncols == 0 || nrows != ncols) {
		printf("Skipping bank %d (%d x %d is absurd or nonsquare)\n",
				current_bank,nrows,ncols);
		return -1;
	}

	/* Get the modulus */

	s	= strstr(t,"Z/");
	if (s==NULL) {
		printf("Skipping bank %d (no modulus found)\n",current_bank);
		return -1;
	}

	s	+= 2;
	stub	 = strcspn(s,"Z ");

	if (*s=='(' && s[stub-1]==')') {
		s++; stub-=2;
	}

	if (stub<=0 || stub>=MODULUSNAME_MAX_LENGTH) {
		printf("Skipping bank %d (bad modulus found)\n",current_bank);
		return -1;
	}

	strncpy(modulusname,s,stub);
	modulusname[stub]=0;

	if (bw_read_modulus_info(modulusname, 1)<0) {
		printf("Skipping bank %d "
			"(modulus information for %s is bad or unavailable)\n",
			current_bank,modulusname);
		return -1;
	}

	/* Everything is allright */

	*p_fm	= fopen(fnm,"r");
	*p_fv	= fopen(fnv,"r");

	if (    *p_fm == NULL || fseek(*p_fm,idlm,SEEK_SET)<0 ||
		*p_fv == NULL || fseek(*p_fv,idlv,SEEK_SET)<0)
	{
		printf("Skipping bank %d (Data disappeared!)\n",current_bank);
		return -1;
		/* We may want to die here */
	}
		
#undef fnm
#undef fnv
#undef fnr

	return 0;
}

