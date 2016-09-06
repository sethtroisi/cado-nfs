#ifndef TREE_STATS_HPP_
#define TREE_STATS_HPP_

#include <string>
#include <vector>
#include <map>

class tree_stats {
    struct function_stats {
        unsigned int ncalled;
        unsigned int min_inputsize;
        unsigned int max_inputsize;
        unsigned int sum_inputsize;
        double spent;
        double projected_time;
        unsigned int projected_calls;
        std::map<std::string, double> small_steps;
        function_stats() : 
            ncalled(0),
            min_inputsize(UINT_MAX),
            max_inputsize(0),
            sum_inputsize(0),
            spent(0),
            small_steps()
        {}
    };

    struct level_stats : public std::map<std::string, function_stats> {
        double projected_time(unsigned int);
        double last_printed_projected_time;
    };


    struct running_stats {
        std::string func;
        unsigned int inputsize;
        double time_self;
        double time_children;
        std::map<std::string, double> small_steps;
        double * substep;
        running_stats() : func()
                          , inputsize(0)
                          , time_self(0)
                          , time_children(0)
                          , small_steps()
                          , substep(NULL)
        { }
    };

    std::vector<level_stats> stack;
    std::vector<running_stats> curstack;

    unsigned int tree_total_breadth;
    double last_print_time;

    double begin;       /* stored as a wct_seconds() return value */


    void print(unsigned int level);

public:
    unsigned int depth;
    
    void enter(const char * func, unsigned int inputsize); 
    int leave(int rc);

    void begin_smallstep(const char * func MAYBE_UNUSED);
    void end_smallstep();
};

#endif	/* TREE_STATS_HPP_ */
