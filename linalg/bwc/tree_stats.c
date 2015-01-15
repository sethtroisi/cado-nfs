#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "portability.h"
#include "select_mpi.h"
#include "utils.h"
#include "tree_stats.h"

void tree_stats_print(tree_stats_ptr stats, unsigned int level)
{
    double sum = 0;
    int nstars=0;
    int nok=0;
    double firstok = 0;
    double complement = 0;
    for(unsigned int k = 1 ; k < 64 ; k++) {
        tree_level_stats_ptr u = stats->stats_stack + k;

        if (k > level && !u->ncalled)
            break;

        if (!u->ncalled) {
            printf("%u *\n", k-1);
            nstars++;
            continue;
        }
        u->projected_time = stats->tree_total_breadth * u->spent / u->sum_inputsize;
        if (nstars && nok == 0) {
            firstok = u->projected_time;
        } else if (nstars && nok == 1) {
            double ratio = firstok / u->projected_time;
            complement = ratio * firstok * (pow(ratio, nstars) - 1) / (ratio - 1);
        }
        nok++;

        sum += u->projected_time;
        printf("%u [%u-%u, %s] %u/%u %.1f -> %.1f (total: %.1f)\n",
                k-1, u->min_inputsize, u->max_inputsize,
                u->func,
                u->ncalled, 1u << (k-1),
                u->spent / u->ncalled, u->projected_time, sum);
        u->last_printed_projected_time = u->projected_time;
    }
    if (nstars && nok >= 2) {
        printf("expected time for levels 0-%u: %.1f (total: %.1f)\n",
                nstars-1, complement, sum + complement);
    }
}

void tree_stats_enter(tree_stats_ptr stats, const char * func, unsigned int inputsize)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { ++stats->depth; return; }
    tree_running_stats_ptr s = stats->stats_curstack + ++stats->depth;
    memset(s, 0, sizeof(struct tree_running_stats_s));
    s->time_self -= wct_seconds();
    s->func = func;
    s->inputsize = inputsize;
}

int tree_stats_leave(tree_stats_ptr stats, int rc)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { --stats->depth; return rc; }
    double now = wct_seconds();
    tree_running_stats_ptr s = stats->stats_curstack + stats->depth;
    s->time_self += now;
    double tt = s->time_self;
    s->time_self -= s->time_children;
    unsigned int level = stats->depth;

    tree_level_stats_ptr t = stats->stats_stack + level;
    if (t->ncalled == 0) {
        t->min_inputsize = s->inputsize;
        t->max_inputsize =s-> inputsize;
        t->func = s->func;
    }
    t->ncalled++;
    if (s->inputsize < t->min_inputsize)
        t->min_inputsize = s->inputsize;
    if (s->inputsize > t->max_inputsize)
        t->max_inputsize = s->inputsize;
    t->sum_inputsize += s->inputsize;
    t->spent += s->time_self;
    if (t->func != s->func && !t->warned) {
        fprintf(stderr, "Warning: at depth level %d (input size %u to %u), function called is not constant: %s or %s\n",
                level, t->min_inputsize, t->max_inputsize, t->func, s->func);
        t->warned=1;
    }

    stats->stats_curstack[stats->depth - 1].time_children += tt;
    --stats->depth;

    /* Is it any useful to print something new ? */

    if (now < stats->last_print_time + 2)
        return rc;

    int needprint = 0;
    for(unsigned int k = 0 ; !needprint && k < 64 ; k++) {
        tree_level_stats_ptr u = stats->stats_stack + k;

        if (k > level && !u->ncalled)
            break;

        if (!u->ncalled)
            continue;

        /* compute the projected time */
        u->projected_time = stats->tree_total_breadth * u->spent / u->sum_inputsize;
        if (u->projected_time < 0.98 * u->last_printed_projected_time)
            needprint = 1;
        else if (u->projected_time > 1.02 * u->last_printed_projected_time)
            needprint = 1;
    }

    if (!needprint)
        return rc;

    stats->last_print_time = now;

    tree_stats_print(stats, level);

    return rc;
}



