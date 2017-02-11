#include "cado.h"
#include "tdict.hpp"
#include "params.h"

int tdict::global_enable = 0;

void tdict_decl_usage(param_list pl)
{
    param_list_decl_usage(pl, "T",   "enable fine-grain timings (use twice to get them for each q)");
}

void tdict_configure_switch(param_list pl)
{
    param_list_configure_switch(pl, "-T", &tdict::global_enable);
}
