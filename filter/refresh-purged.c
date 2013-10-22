
/* This tool may be used to transform an old style .nodup + .purged
 * combination into a new style .purged file. Use as follows:
 *
 *  ./refresh-purged -oldpurged /localdisk/tmp/160/160.purged  -oldnodup /localdisk/tmp/160/160.nodup.gz  -out /localdisk/tmp/160/160.purged2
 *
 */

#include "cado.h"
#include <stdio.h>
#include <string.h>
#include "portability.h"
#include "utils.h"

int
main (int argc, char **argv)
{
    param_list pl;
    param_list_init(pl);

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        abort(); 
    }

    const char * oldp = param_list_lookup_string(pl, "oldpurged");
    const char * oldn = param_list_lookup_string(pl, "oldnodup");
    const char * outname = param_list_lookup_string(pl, "out");

    ASSERT_ALWAYS(oldp && oldn && outname);

    if (param_list_warn_unused(pl))
        exit(1);

    /* Using a purgedfile_stream is going to work, even though the
     * recognized (a,b) values will be complete rubbish.
     */
    char line[RELATION_MAX_BYTES];
    purgedfile_stream ps;
    purgedfile_stream_init(ps);
    purgedfile_stream_openfile(ps, oldp);
    ps->parse_only_ab = 1;

    relation_stream rs;
    relation_stream_init(rs);
    relation_stream_openfile(rs, oldn);
    rs->parse_only_ab = 1;

    FILE * out = fopen(outname, "w");
    fprintf(out, "%" SCNu64 " %" SCNu64 "\n", ps->nrows, ps->ncols);

    int rn = 0;
    for( ; purgedfile_stream_get(ps, line) >= 0 ; ) {
        ASSERT_ALWAYS(rn <= ps->nodup_index);
        for( ; rn < ps->nodup_index ; rn++) {
            ASSERT_ALWAYS(relation_stream_get(rs, NULL, 0, 0, 10, 0));
            if (relation_stream_disp_progress_now_p(rs)) {
                fprintf(stderr, "Read %" PRIu64 " relations in %.1f s -- "
                                "%.1f MB/s\n", rs->nrels, rs->dt, rs->mb_s);
            }
        }
        char * p;
        ASSERT_ALWAYS((p = strchr(line, ' ')) != NULL);
        fprintf(out,"%d %" PRId64 " %" PRIu64 "%s", ps->nodup_index, rs->rel.a, rs->rel.b, p);
        if (purgedfile_stream_disp_progress_now_p(ps)) {
            fprintf(stderr, "Treated %d/%" PRIu64 " purged relations\n",
                    ps->rrows, ps->nrows);
        }
    }
    fclose(out);
    relation_stream_closefile(rs);
    relation_stream_clear(rs);
    purgedfile_stream_closefile(ps);
    purgedfile_stream_clear(ps);
    return 0;
}
