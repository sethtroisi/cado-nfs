#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "filenames.h"
#include "bw_scalar.h"
#include "variables.h"
#include "tagfile.h"

int read_tag_file(void)
{
	FILE *f;
	int res;
	char modulusname[MODULUSNAME_MAX_LENGTH];
	f=fopen(wip_tag_filename,"r");
	if (f == NULL)
		die("%s : %s\n",1,wip_tag_filename,strerror(errno));
	res=(fscanf(f,":%d:%d:%d:%d:%d:" MODULUSNAME_FMT "\n",
				&computed_m_param,&computed_n_param,&ncols,
				&total_work,&periodicity,modulusname)==6);
	nrows=ncols;
	fclose(f);
	if (!res)
		die("Unable to read info from %s\n",1,wip_tag_filename);

	if (bw_read_modulus_info(modulusname, 1)<0)
		die("Unable to read info for modulus %s\n",1,modulusname);

	return 0;
}

int write_tag_file(char * modulusname)
{
	FILE	      * f;
	int		s,t;

	s = t = ncols+2*n_param;

	s = (s / n_param) + ( (s % n_param) != 0);
	t = (t / m_param) + ( (t % m_param) != 0);

	total_work = s + t;

	if (periodicity==0)
		periodicity = MAX(100 , total_work / 500);
	
	f = fopen(wip_tag_filename,"w");
	if (f == NULL)
		die("%s : %s",1,wip_tag_filename,strerror(errno));
	
	fprintf(f,":%d:%d:%d:%d:%d:%s\n",
			m_param,n_param,ncols,
			total_work,periodicity,modulusname);
	
	fprintf(f,"# DO NOT DELETE THIS FILE. The `slaves' rely on it\n");
	fclose(f);

	return 0;
}

void load_x_vectors(coord_t * tab, int begin, int number)
{
	int i;
	char name[FILENAME_LENGTH];
	FILE *f;

	for(i=0;i<number;i++) {
		sprintf(name,x_meta_filename, begin + i);
		f=fopen(name,"r");
		if (f==NULL)
			die("fopen(%s) : %s\n",1,name,strerror(errno));
		if (fscanf(f,"e%d\n",&tab[i])!=1)
			die("%s : no data\n",1,name);
		fclose(f);
	}
}
