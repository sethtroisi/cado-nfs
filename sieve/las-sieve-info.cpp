#include "cado.h"

#include "las-sieve-info.hpp"

/* This function creates a new sieve_info structure, taking advantage of
 * structures which might already exist in las.sievers
 *  - for sieving, if factor base parameters are similar, we're going to
 *    share_factor_bases()
 *  - for cofactoring, if large prime bounds and mfbs are similar, we're
 *    going to reuse the strategies.
 *
 * The siever_config structure to be passed to this function is not
 * permitted to lack anything.
 *
 * This function differs from las_info::get_sieve_info_from_config(),
 * since the latters also registers the returned object within
 * las.sievers (while the function here only *reads* this list).
 */
sieve_info::sieve_info(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & sievers, cxx_param_list & pl, bool try_fbc) /*{{{*/
    : cpoly(cpoly), conf(sc)
{
    I = 1UL << sc.logI_adjusted;

    std::list<sieve_info>::iterator psi;

    /*** Sieving ***/

    psi = find_if(sievers.begin(), sievers.end(), sc.same_fb_parameters());

    if (psi != sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy factor base data from previous siever\n");
        sides[0].fb = other.sides[0].fb;
        sides[1].fb = other.sides[1].fb;
    } else {
        const char * fbc_filename = param_list_lookup_string(pl, "fbc");
        if (try_fbc && fbc_filename) {
            FILE * f = fopen(fbc_filename, "r");
            sides[0].fb=std::make_shared<fb_factorbase>(f, sc, 0);
            sides[1].fb=std::make_shared<fb_factorbase>(f, sc, 1);
            /* what shall we do ? If we copy, instead of mmaping, then
             * we'll duplicate the memory, so it's no good. OTOH, this
             * stuff now has many many pointers...
             */
            fclose(f);
        } else {
            /* TODO: get rid of the indirection. We can't have sieve_info
             * contain a full-fledged reference because then the type
             * would not be default constructible.
             */
            sides[0].fb=std::make_shared<fb_factorbase>(*(cxx_cado_poly const *)cpoly, sc, pl, 0);
            sides[1].fb=std::make_shared<fb_factorbase>(*(cxx_cado_poly const *)cpoly, sc, pl, 1);
        }
        /* we used to call print_fb_statistics at this point, it should
         * no longer be needed now.  */
    }


    /*** Cofactoring ***/
    psi = find_if(sievers.begin(), sievers.end(), sc.same_cofactoring());

    if (psi != sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy cofactoring strategies from previous siever\n");
        strategies = other.strategies;
    } else {
        init_strategies(pl);
    }
}
/*}}}*/

void sieve_info::update_norm_data ()/*{{{*/
{
    for(int side = 0 ; side < 2 ; side++) {
        sides[side].lognorms = std::make_shared<lognorm_smart>(conf, cpoly, side, qbasis, J);
    }
}

/*}}}*/
void sieve_info::update (unsigned int nr_workspaces)/*{{{*/
{
    uint64_t A = UINT64_C(1) << conf.logA;

    for(int side = 0 ; side < 2 ; side++) {
        sieve_info::side_info & sis(sides[side]);
        if (!sis.fb) continue;

        /* TODO: now we can adjust this according to logI !!! 
         *
         * That's actually the whole point of this overhaul !
         */
        fbprime_t bk_thresh = conf.bucket_thresh;
        fbprime_t bk_thresh1 = conf.bucket_thresh1;
        fbprime_t fbb = conf.sides[side].lim;
        if (bk_thresh > fbb) bk_thresh = fbb;
        if (bk_thresh1 == 0 || bk_thresh1 > fbb) bk_thresh1 = fbb;

        fb_factorbase::key_type K {
            {bk_thresh, bk_thresh1, fbb, fbb},
            conf.td_thresh,
            conf.skipped,
            sis.lognorms->scale,
            nr_workspaces
        };

        /* The size of the slices must be in accordance to our
         * multithread setting. Hopefully, that one does not change...
         * We want to divide in small enough pieces, so that the amout of
         * work sums up more or less evenly across the threads.
         *
         * This "small enough" criterion used to be computed around this
         * place, and has now moved to inside the make_slices call. It is
         * inherently tied to the count of the entries in each part.
         */

        sis.fbs = &(*sis.fb)[K];
    }

    /* TODO: toplevel depends on the fb thresholds. Or more specifically,
     * it is decided once we have interpreted the fb thresholds, and
     * decided on where are the different parts.
     *
     * Now since the fb_thresholds are interpreted at the same time as
     * the log scale parameter, we must be done with the log scale
     * computation before we do anything top-level related.
     */

    // Now that fb have been initialized, we can set the toplevel.
    toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        if (!sides[side].fb) continue;
        int level = sides[side].fbs->get_toplevel();
        if (level > toplevel) toplevel = level;
    }

#if 0
    /* Initialize the number of buckets */
    /* (it's now done in sieve_info::update, which is more timely) */
    /* set the maximal value of the number of buckets. This directly
     * depends on A */
    uint32_t XX[FB_MAX_PARTS] = { 0, NB_BUCKETS_2, NB_BUCKETS_3, 0};
    uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
    uint64_t A = UINT64_C(1) << conf.logA;
    XX[toplevel] = iceildiv(A, BRS[toplevel]);
    for (int i = toplevel+1; i < FB_MAX_PARTS; ++i)
        XX[i] = 0;
    if (toplevel > 1 && XX[toplevel] == 1) {
        XX[toplevel-1] = iceildiv(A, BRS[toplevel-1]);
        ASSERT_ALWAYS(XX[toplevel-1] != 1);
    }
    for (int i = 0; i < FB_MAX_PARTS; ++i) {
        nb_buckets[i] = XX[i];
    }
#endif

    for(int side = 0 ; side < 2 ; side++) {
	init_trialdiv(side); /* Init refactoring stuff */
        if (!sides[side].fb) continue;
        init_fb_smallsieved(side);
        /*
        verbose_output_print(0, 2, "# small side-%d factor base", side);
        size_t nr_roots;
        sides[side].fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
        */
    }

    /* update number of buckets at toplevel */
    size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
    nb_buckets[toplevel] = iceildiv(A, BRS[toplevel]);

    // maybe there is only 1 bucket at toplevel and less than 256 at
    // toplevel-1, due to a tiny J.
    if (toplevel > 1) {
        if (nb_buckets[toplevel] == 1) {
            nb_buckets[toplevel-1] = iceildiv(A, BRS[toplevel - 1]);
            // we forbid skipping two levels.
            ASSERT_ALWAYS(nb_buckets[toplevel-1] != 1);
        } else {
            nb_buckets[toplevel-1] = BRS[toplevel]/BRS[toplevel-1];
        }
    }
}/*}}}*/


void sieve_info::share_factor_bases(sieve_info& other)
{
    for(int side = 0 ; side < 2 ; side++) {
        ASSERT_ALWAYS(conf.sides[side].lim == other.conf.sides[side].lim);
        ASSERT_ALWAYS(conf.sides[side].powlim == other.conf.sides[side].powlim);
        sieve_info::side_info & sis(sides[side]);
        sieve_info::side_info & sis0(other.sides[side]);
        sis.fb = sis0.fb;
    }
}
