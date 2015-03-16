#include "cado.h"
#include "las-dlog-base.h"
#include "las-types.h"

#include <cstdio>
#include <cstdint>
#include <cctype>


#include "utils.h"

void las_dlog_base::declare_parameter_usage(param_list pl)
{
    param_list_decl_usage(pl, "renumber", "renumber table (for the descent)");
    param_list_decl_usage(pl, "log", "log table, as built by reconstructlog");
}

void las_dlog_base::lookup_parameters(param_list pl)
{
    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "renumber")) != NULL) {
        renumberfilename = strdup(tmp);
    }
    if ((tmp = param_list_lookup_string(pl, "log")) != NULL) {
        logfilename = strdup(tmp);
    }
    if (!logfilename != !renumberfilename) {
        fprintf(stderr, "In descent mode, want either renumber+log, or none\n");
        exit(EXIT_FAILURE);
    }
}

las_dlog_base::las_dlog_base()
{
    renumberfilename = NULL;
    logfilename = NULL;
    renumber_init_for_reading(renumber_table);
}

las_dlog_base::~las_dlog_base()
{
    renumber_clear (renumber_table);
    free(renumberfilename);
    free(logfilename);
}

void las_dlog_base::set_default_lpb(siever_config_srcptr sc)
{
    lpb[0] = sc->sides[0]->lpb;
    lpb[1] = sc->sides[1]->lpb;
}

bool las_dlog_base::is_known(int side, uint64_t p, uint64_t r) const {
    if (renumberfilename) {
        index_t h = renumber_get_index_from_p_r(renumber_table, p, r, side);
        return known_logs[h];
    } else {
        return p <= (1UL<<lpb[side]);
    }
}


void las_dlog_base::read()
{
    if (!renumberfilename) {
        verbose_output_print(0, 1, "# Descent: no access to renumber table given, using lpb(%lu/%lu) to decide what are the supposedly known logs\n",
                lpb[0], lpb[1]);
        return;
    }

    verbose_output_print(0, 1, "# Descent: will get list of known logs from %s, using also %s for mapping\n", logfilename, renumberfilename);

    renumber_read_table(renumber_table, renumberfilename);
    uint64_t nprimes = renumber_table->size;
    known_logs.reserve(nprimes);
    /* format of the log table: there are FIVE different line types.
     *
     * [index] added column [log]
     * [index] bad ideals [log]
     * [index] [p] [side] rat [log]
     * [index] [p] [side] [r] [log]
     * [index] SM col [i] [log]
     *
     * Here we care only about the index anyway. By the way, this index
     * is written in hex.
     */
    FILE * f = fopen(logfilename, "r");
    ASSERT_ALWAYS(f != NULL);
    for(int lnum = 0 ; ; lnum++) {
        char line[1024];
        char * x = fgets(line, sizeof(line), f);
        if (x == NULL) break;
        for( ; *x && isspace(*x) ; x++) ;
        if (*x == '#') continue;
        if (!*x) continue;
        errno=0;
        unsigned long z = strtoul(x, &x, 16);
        if (errno) {
            fprintf(stderr, "Parse error at line %d in %s: %s\n", lnum, logfilename, line);
            break;
        }
        known_logs[z] = true;
    }
    fclose(f);
}


