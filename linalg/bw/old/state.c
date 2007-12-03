#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif

#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "macros.h"
#include "auxfuncs.h"
#include "types.h"
#include "macros.h"
#include "filenames.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "options.h"
#include "corners.h"
#include "slave.h"
#include "variables.h"
#include "bw_lvblock_steps.h"
#include "threaded.h"
#include "barrier.h"
/* #include "version.h" */
#include "certificate.h"
#include "recovery.h"

static struct vs_data ** vsd_table;

/* Huge warning ! There must be a synchronization step before this
 * function.
 */
void commit_more_values(bw_vector_block w, bw_scalar t)
{
	int i,l;
	int d,s;
	s=0;
	for(i=tsv()->x0;i<tsv()->x1;i++) {
		d=bw_x[i]-s;
		s=bw_x[i];
		bw_lvblock_step_n00(w,d);
		for(l=0;l<nbys;l++) {
			/* NO ! Don't do this, since threads might be
			 * currently be reading at this location (it's the
			 * source vector of the following operation). The
			 * reduction has already been performed anyway in
			 * lvblock_reduce */
			/* bw_reduce_short_scalar(w); */
			bw_scalar_write(vsd_table[l]->a_save[i],w);
			bw_lvblock_step_00p(w);
		}
		bw_lvblock_step_00z(w);
	}
}


void state_checkpoint(bw_vector_block v, int force_certif)
{
	int	j;

	for(j=0;j<nbys;j++) {
		vsd_table[j]->vec=v;
		vs_checkpoint(vsd_table[j], force_certif);
	}
}

int recover_state(bw_vector_block v, struct mksol_info_block * msi)
{
	int	j;

	vsd_table=malloc(nbys*sizeof(struct vs_data *));
	for(j=0;j<nbys;j++) {
		vs_try_startup_sync(&(vsd_table[j]),first.j+j,
				msi->differential);
	}

	for(j=1;j<nbys;j++) {
		if (vsd_table[j]->rev_v != vsd_table[0]->rev_v) {
			die("Inconsistent revisions reached\n",1);
		}
	}

	for(j=0;j<nbys;j++) {
		vsd_table[j]->msi=msi;
		vsd_table[j]->sum=msi->sum;
		vsd_table[j]->vec=v;
		vs_do_startup_sync(vsd_table[j]);
		bw_lvblock_step_00p(v);
	}

	printf("Recovered: %d * %d = %d iterations\n",
			vsd_table[0]->rev_v,periodicity,
			vsd_table[0]->rev_v * periodicity);

	return vsd_table[0]->rev_v * periodicity;
}

/* vim:set sw=8: */
