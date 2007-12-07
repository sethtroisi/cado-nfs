#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "old-endian.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "variables.h"
#include "options.h"
#include "filenames.h"
#include "tagfile.h"

void showuse(void)
{
	die("Usage : bw-addvec <bank#> [<vectors> ...]\n",1);
}

int main(int argc, char * argv[])
{
	FILE ** fds;
	bw_scalar s,*x;
	int i;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	if (argc <= 2) {
		showuse();
	}

	check_endianness(NULL);
	bank_num=atoi(argv[1]);
	set_all_filenames(bank_num);
	read_tag_file();


	if (isatty(fileno(stdout))) {
		fprintf(stderr,"Won't dump to a tty\n");
		showuse();
	}

	fds = malloc((argc-2)*sizeof(FILE*));
	for(i=0;i<argc-2;i++) {
		fds[i]=fopen(argv[i+2],"r");
		if (fds[i]==NULL) {
			perror(argv[i+2]);
			showuse();
		}
	}
	bw_scalar_alloc(s,BW_SCALAR_LONG);
	x=malloc((argc-2)*sizeof(bw_scalar));
	for(i=0;i<argc-2;i++) {
		bw_scalar_alloc(x[i],BW_SCALAR_LONG);
	}

	for(;;) {
		int nread=0;
		memset(s,0,bw_longsize*sizeof(mp_limb_t));
		for(i=0;i<argc-2;i++) {
			memset(x[i],0,bw_longsize*sizeof(mp_limb_t));
			nread+=bw_scalar_read(x[i],fds[i]);
		}
		if (nread < (argc-2)*bw_filesize) {
			if (nread) {
				fprintf(stderr,"Wrong size\n");
			}
			break;
		}
		bw_long_scalar_set_zero(s);
		for(i=0;i<argc-2;i++) {
			mpn_add_n(s,s,x[i],bw_longsize);
		}
		mpn_tdiv_qr(s+bw_allocsize+2,s,
				0,
				s, bw_allocsize+2,
				modulus_plain, bw_allocsize);
		bw_scalar_write(stdout,s);
	}

	for(i=0;i<argc-2;i++) {
		free(x[i]);
		fclose(fds[i]);
	}
	return 0;
}
/* vim:set sw=8: */
