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
#include "utils.h"
#include "filter_utils.h"
#include "gzip.h"
#include "filter_matrix.h"
#include "sparse.h"
#include "report.h"

/* terrific hack: everybody on the same line
   the first is to be destroyed in replay!!! */
void
reportn (report_t *rep, int32_t *ind, int n, MAYBE_UNUSED int32_t j)
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
    if (j >= 0)
        fprintf(rep->outfile, " #%d", j);
#endif
	fprintf(rep->outfile, "\n");
    }
    else if((rep->type == 1) && (rep->mark != -2)){
	// mark == -2 => we are probably resuming and we don't care
	// to consume a lot of memory that will not be used, anyway.
	rep->mark += 1;
	if(rep->history[rep->mark] == NULL)
	    rep->history[rep->mark] = 
		(int32_t *)malloc((rep->bufsize+1)*sizeof(int32_t));
	rep->history[rep->mark][0] = n;
	for(i = 0; i < n; i++)
	    rep->history[rep->mark][i+1] = ind[i];
    }
}

/* print a new line "i" in the history file */
void
report1 (report_t *rep, int32_t i, int32_t j)
{
  reportn (rep, &i, 1, j);
}

/* print a new line "i1 i2" in the history file */
void
report2 (report_t *rep, int32_t i1, int32_t i2, int32_t j)
{
    int32_t tmp[2];

    tmp[0] = i1;
    tmp[1] = i2;
    reportn (rep, tmp, 2, j);
}
