#include "cado.h"
#include <pthread.h>
#include "tdict.hpp"
#include "params.h"

namespace tdict {

int global_enable = 0;

#ifndef DISABLE_TIMINGS

pthread_mutex_t slot_base::m = PTHREAD_MUTEX_INITIALIZER;

void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "T",   "enable fine-grain timings (use twice to get them for each q)");
}

void configure_aliases(cxx_param_list &)
{
}

void configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-T", &global_enable);
}

#else

void declare_usage(cxx_param_list &) {}
void configure_aliases(cxx_param_list &) {}
void configure_switches(cxx_param_list &) {}

#endif
};
