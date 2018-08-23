#ifndef LAS_PARALLEL_HPP_
#define LAS_PARALLEL_HPP_

#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif
#include <exception>
#include <memory>
#include "params.h"

class las_parallel_desc {
    const char * desc_c = NULL;
    double jobram = -1;

    int subjobs_per_cpu_binding_zone = 1;
    int threads_per_subjob = 1;
    int memory_binding_zones = 1;
    int memory_binding_size = 0;       /* unset in the default setting */
    int cpu_binding_zones_per_memory_binding_zone = 1;
    int cpu_binding_size = 0;

    /* when we apply the "loose" cpu binding, we often do this to process
     * or prepare data that is of interest to all cores, and we wish to
     * do that collectively. We wonder, however, how many threads this
     * means. The following rules are followed:
     *
     *  - no hwloc: same integer as the number of threads of the unique
     *  subjob.
     *  - without job replication: total number of PUs in the binding context.
     *  - with job replication: total number of PUs on the machine.
     */
    int number_of_threads_loose() const;

    struct helper;
    friend struct helper;
    std::shared_ptr<helper> help;
public:
    int number_of_memory_binding_zones() const { return memory_binding_zones; }
    int number_of_subjobs_per_cpu_binding_zone() const { return subjobs_per_cpu_binding_zone; }
    int number_of_threads_per_subjob() const { return threads_per_subjob; }

    int number_of_subjobs_per_memory_binding_zone() const {
        return subjobs_per_cpu_binding_zone *
            cpu_binding_zones_per_memory_binding_zone;
    }
    int number_of_subjobs_total() const {
        return memory_binding_zones *
            number_of_subjobs_per_memory_binding_zone();
    }
    int number_of_threads_total() const {
        return number_of_subjobs_total() * threads_per_subjob;
    }
    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "t", "Number of threads and subjobs. Use -t help for extended documentation");
#ifdef HAVE_HWLOC
        param_list_decl_usage(pl, "job-memory", "Estimated memory per subjobs, used for job placement (see -t help)");
#else
        param_list_decl_usage(pl, "job-memory", "(unused, needs hwloc)\n");
#endif
    }
    las_parallel_desc() = default;
    las_parallel_desc(las_parallel_desc const &) = default;
    las_parallel_desc(cxx_param_list & pl, double jobram = -1);
    void display_binding_info() const ;
    int set_loose_binding() const ;
    int set_subjob_binding(int k) const ;

    struct needs_job_ram : public std::exception {
        virtual const char * what() const noexcept {
            return "The \"fit\" specifier requires an estimate of the available RAM\n";
        }
    };

    class bad_specification : public std::exception {
        std::string message;
        void build_message(std::ostream &) { }
        template <typename Car, typename... Cdr>
        void build_message(std::ostream & os, Car car, Cdr... cdr) {
            build_message(os << car, cdr...);
        }
        public:
        template<typename ...Args>
        bad_specification(Args&&... args) {
            std::ostringstream os;
            build_message(os, args...);
            message = os.str();
        }
        virtual const char * what() const noexcept { return message.c_str(); }
    };
};

#endif	/* LAS_PARALLEL_HPP_ */
