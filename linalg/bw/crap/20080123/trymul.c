#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "old-endian.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "bw_lvblock_steps.h"
#include "variables.h"
#include "slave.h"
#include "filenames.h"
#include "matrix.h"
#include "tagfile.h"
#include "timer.h"
#include "threaded.h"
#include "addmul.h"
/* #include "version.h" */
#include "gmp-hacks.h"

int computed_nbxs=1;
int computed_nbys=1;

coord_t *bw_x;		/* These are just elements of the canonical basis */

void showuse(void)
{
	die("Usage : bw-trymul <bank#>\n",1);
}

int main(int argc, char *argv[])
{
	bw_vector_block bw_v, bw_w;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	/* fputs(version_string,stderr); */
	check_endianness(stderr);

	if (argc!=2)
		showuse();
	bank_num=atoi(argv[1]);

	set_all_filenames(bank_num);
	read_tag_file();

#ifdef HARDCODE_PARAMS
	consistency_check("nbys",computed_nbys,nbys);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	load_matrix();

	configure_threads(1,ncols);
	compute_thread_offsets();

	bw_lvblock_alloc(bw_v);
	bw_lvblock_set_zero_separated(bw_v);
	bw_lvblock_read(bw_v,stdin);
	bw_lvblock_reduce_separated(bw_v);
	bw_lvblock_alloc(bw_w);
	bw_lvblock_set_zero_separated(bw_w);

	multiply(bw_w,bw_v);
	bw_lvblock_reduce_separated(bw_w);
	bw_lvblock_write(stdout,bw_w);

	return 0;
}

/* vim:set sw=8: */
