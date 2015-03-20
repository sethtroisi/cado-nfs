#ifndef LAS_DLOG_BASE_H_
#define LAS_DLOG_BASE_H_

#include "utils.h"
#include "las-types.h"

#ifndef __cplusplus
#error "This file is C++ only"
#endif

#include <vector>

struct las_dlog_base {
private:
    char * renumberfilename;
    char * logfilename;

    renumber_t renumber_table;
    std::vector<bool> known_logs;
    unsigned long lpb[2];

public:
    bool is_known(int side, uint64_t p, uint64_t r) const;
    void lookup_parameters(param_list pl);
    las_dlog_base();
    ~las_dlog_base();
    void read();
    void set_default_lpb(siever_config_srcptr sc);

    static void declare_parameter_usage(param_list pl);
};

#endif	/* LAS_DLOG_BASE_H_ */
