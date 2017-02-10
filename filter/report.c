/* report.c --- auxiliary merge program

Copyright 2008-2009 Francois Morain, Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "filter_config.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "report.h"

/* terrific hack: everybody on the same line
   the first is to be destroyed in replay!!! */
void
reportn (report_t *rep, index_signed_t *ind, int n, MAYBE_UNUSED index_t j)
{
    int i;

#if DEBUG >= 1
    fprintf(stderr, "Reporting for n=%d\n", n);
#endif
    if(rep->type == 2)
	return;
    if(rep->type == 0){
	for(i = 0; i < n; i++){
	    fprintf(rep->outfile, "%ld", (long int) ind[i]);
	    if(i < n-1)
		fprintf(rep->outfile, " ");
	}
#ifdef FOR_DL
        fprintf(rep->outfile, " #%lu", (unsigned long) j);
#endif
	fprintf(rep->outfile, "\n");
    }
    else if((rep->type == 1) && (rep->mark != -2)){
	// mark == -2 => we are probably resuming and we don't care
	// to consume a lot of memory that will not be used, anyway.
	rep->mark += 1;
	if(rep->history[rep->mark] == NULL)
	    rep->history[rep->mark] = 
		(index_t *)malloc((rep->bufsize+1)*sizeof(index_t));
	rep->history[rep->mark][0] = n;
	for(i = 0; i < n; i++)
	    rep->history[rep->mark][i+1] = ind[i];
    }
}

/* Print a new line "i" in the history file.
   This function is always called with i >= 0. */
void
report1 (report_t *rep, index_signed_t i, index_t j)
{
  ASSERT(i >= 0);
  reportn (rep, &i, 1, j);
}
