#include "cado.h"

#include "las-sieve-shared-data.hpp"
#include "las-cofactor.hpp"     // proxy facul_make_strategies


void sieve_shared_data::declare_usage(cxx_param_list & pl)
{
    cxx_cado_poly::declare_usage(pl);
    param_list_decl_usage(pl, "fb0",   "factor base file on the rational side");
    param_list_decl_usage(pl, "fb1",   "factor base file on the algebraic side");
    param_list_decl_usage(pl, "fbc",  "factor base cache file (not yet functional)");
}

sieve_shared_data::side_data::side_data(int side,
        cxx_cado_poly const & cpoly,
        cxx_param_list & pl)
    :
        f(cpoly->pols[side]),
        fb(cpoly, side, pl, param_list_lookup_string(pl, "fbc"))
{
}

/* FIXME: get_trialdiv_data is currently in las-trialdiv.cpp ; it belongs
 * here (that is, the bulk should stay there, but the interaction with
 * the cache mechanism should be be here instead).
 */

/* The fb_factorbase ctor parses the lim[01] and powlim[01] directly from
 * the command line. These get interpreted as the "max" bounds, and we
 * compute the factor base up to that limit.
 */
sieve_shared_data::sieve_shared_data(cxx_cado_poly const & cpoly, cxx_param_list & pl) /*{{{*/
    : cpoly(cpoly), sides {
        {0, cpoly, pl},
        {1, cpoly, pl}
    }
{
    cofactfilename = param_list_lookup_string (pl, "file-cofact");
}
/*}}}*/
unsieve_data const * sieve_shared_data::get_unsieve_data(siever_config const & conf) /* {{{ */
{
    std::pair<int, int> p(conf.logI, conf.logA);
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    auto it = us_cache.find(p);
    if (it != us_cache.end()) {
        pthread_mutex_unlock(&lock);
        return &it->second;
    }
    auto itb = us_cache.insert(std::make_pair(p, unsieve_data(p)));
    ASSERT(itb.second);
    pthread_mutex_unlock(&lock);
    return &(*itb.first).second;
}/*}}}*/
j_divisibility_helper const * sieve_shared_data::get_j_divisibility_helper(int J) /* {{{ */
{
    ASSERT_ALWAYS(J);
    /* Round to the next power of two */
    unsigned int Jround = 1 << nbits(J-1);
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    auto it = jdiv_cache.find(Jround);
    if (it != jdiv_cache.end()) {
        pthread_mutex_unlock(&lock);
        return &it->second;
    }
    auto itb = jdiv_cache.insert(std::pair<unsigned int, j_divisibility_helper>(Jround, Jround));
    ASSERT(itb.second);
    pthread_mutex_unlock(&lock);
    return &(*itb.first).second;
}/*}}}*/
facul_strategies_t const * sieve_shared_data::get_strategies(siever_config const & conf) /* {{{ */
{
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    auto it = facul_strategies_cache.find(conf);
    if (it != facul_strategies_cache.end()) {
        pthread_mutex_unlock(&lock);
        return it->second.get();
    }


#if 0
    /* Temporary hack. We return *ALWAYS* the same cofactoring strategy.
     * TODO: investigate, see how this behaves. (for the descent case)
     */
    /* Cannot work: the descent is _really_ allowed to have various mfb
     * set up, so that multiple strategy tables are necessary.
     */
    if (!facul_strategies_cache.empty()) {
        it = facul_strategies_cache.begin();
        if (!siever_config::has_same_cofactoring(conf)(it->first)) {
            verbose_output_print(0, 1, "# NOTE: using previously stored cofactoring strategy, although it was not necessarily meant for the current set of parameters.\n");
        }
        return it->second.get();
    }
#endif

    FILE* file = NULL;
    if (cofactfilename != NULL) /* a file was given */
        file = fopen (cofactfilename, "r");
    double time_strat = seconds();

    auto itb = facul_strategies_cache.insert(std::make_pair(conf,
            std::shared_ptr<facul_strategies_t>(
                facul_make_strategies (conf, file, 0),
                facul_clear_strategies))
            );

    ASSERT_ALWAYS(itb.second);
    verbose_output_print(0, 1, "# Building/reading strategies took %1.1fs\n",
            seconds() - time_strat);

    if (file)
        fclose (file);

    if (!(*itb.first).second.get()) {
        fprintf (stderr, "impossible to read %s\n", cofactfilename);
        pthread_mutex_unlock(&lock);
        abort ();
    }
    pthread_mutex_unlock(&lock);

    return (*itb.first).second.get();
}/*}}}*/

sieve_shared_data::~sieve_shared_data()
{
    char buf1[16];
    verbose_output_print(0, 2, "# Getting rid of sieve_shared_data structure [rss=%s]\n",
            size_disp_fine(1024UL * Memusage2(), buf1, 10000.0));
}
