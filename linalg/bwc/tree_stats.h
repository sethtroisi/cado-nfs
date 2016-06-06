#ifndef TREE_STATS_H_
#define TREE_STATS_H_

/* statistics for projected running time. */
struct tree_level_stats_s {
    unsigned int ncalled;
    unsigned int min_inputsize;
    unsigned int max_inputsize;
    unsigned int sum_inputsize;
    double spent;
    double projected_time;;
    const char * func;
    int warned;
    double last_printed_projected_time;
};
typedef struct tree_level_stats_s tree_level_stats[1];
typedef struct tree_level_stats_s *tree_level_stats_ptr;

struct tree_running_stats_s {
    const char * func;
    unsigned int inputsize;
    double time_self;
    double time_children;
};

typedef struct tree_running_stats_s tree_running_stats[1];
typedef struct tree_running_stats_s *tree_running_stats_ptr;



struct tree_stats_s {
    struct tree_level_stats_s stats_stack[64];
    struct tree_running_stats_s stats_curstack[64];
    unsigned int depth;

    unsigned tree_total_breadth;
    double last_print_time;
};
typedef struct tree_stats_s tree_stats[1];
typedef struct tree_stats_s *tree_stats_ptr;
typedef const struct tree_stats_s *tree_stats_srcptr;


#ifdef __cplusplus
extern "C" {
#endif

void tree_stats_print(tree_stats_ptr stats, unsigned int level);

void tree_stats_enter(tree_stats_ptr stats, const char * func, unsigned int inputsize);

int tree_stats_leave(tree_stats_ptr stats, int rc);

#ifdef __cplusplus
}
#endif

#endif	/* TREE_STATS_H_ */
