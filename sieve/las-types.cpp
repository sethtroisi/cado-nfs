#include "cado.h"
#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>
#include "las-types.hpp"
#include "las-config.h"
#include "las-norms.hpp"

/*****************************/

/* sieve_info stuff */

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
sieve_info::sieve_info(las_info & las, siever_config const & sc, param_list pl)/*{{{*/
{
    cpoly = las.cpoly;
    conf = sc;

    I = 1UL << sc.logI_adjusted;

    std::list<sieve_info>::iterator psi;
    
    /*** Sieving ***/

    psi = find_if(las.sievers.begin(), las.sievers.end(), sc.same_fb_parameters());

    if (psi != las.sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy factor base data from previous siever\n");
        share_factor_bases(other);
    } else {
        verbose_output_print(0, 1, "# bucket_region = %" PRIu64 "\n",
                BUCKET_REGION);
        init_factor_bases(pl);
        for (int side = 0; side < 2; side++) {
            print_fb_statistics(side);
        }
    }

    // Now that fb have been initialized, we can set the toplevel.
    toplevel = MAX(sides[0].fb->get_toplevel(), sides[1].fb->get_toplevel());

    /* If LOG_BUCKET_REGION == sc.logI, then one bucket (whose size is the
     * L1 cache size) is actually one line. This changes some assumptions
     * in sieve_small_bucket_region and resieve_small_bucket_region, where
     * we want to differentiate on the parity on j.
     */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= conf.logI_adjusted);

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
        init_fb_smallsieved(side);
        verbose_output_print(0, 2, "# small side-%d factor base", side);
        size_t nr_roots;
        sides[side].fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
    }

    /*** Cofactoring ***/

    psi = find_if(las.sievers.begin(), las.sievers.end(), sc.same_cofactoring());

    if (psi != las.sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy cofactoring strategies from previous siever\n");
        strategies = other.strategies;
    } else {
        init_strategies(pl);
    }
}
/*}}}*/

void sieve_info::update (size_t nr_workspaces)/*{{{*/
{
#ifdef SMART_NORM
        /* Compute the roots of the polynomial F(i,1) and the roots of its
         * inflection points d^2(F(i,1))/d(i)^2. Used in
         * init_smart_degree_X_norms_bucket_region_internal.  */
        for(int side = 0 ; side < 2 ; side++) {
            if (cpoly->pols[side]->deg >= 2)
                init_norms_roots (*this, side);
        }
#endif

  /* update number of buckets at toplevel */
  uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;

  /* wondering whether having the "local A" at hand would be a plus. */
  uint64_t A = ((uint64_t)J) << conf.logI_adjusted;

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

  /* Update the slices of the factor base according to new log base */
  for(int side = 0 ; side < 2 ; side++) {
      /* The safety factor controls by how much a single slice should fill a
         bucket array at most, i.e., with .5, a single slice should never fill
         a bucket array more than half-way. */
      const double safety_factor = .5;
      sieve_info::side_info & sis(sides[side]);
      double max_weight[FB_MAX_PARTS];
      for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++) {
        max_weight[i_part] = sis.max_bucket_fill_ratio[i_part] / nr_workspaces
            * safety_factor;
      }
      sis.fb->make_slices(sis.scale * LOG_SCALE, max_weight);
  }
}/*}}}*/

/* las_info stuff */

#ifdef  DLP_DESCENT
void las_info::init_hint_table(param_list_ptr pl)/*{{{*/
{
    const char * filename = param_list_lookup_string(pl, "descent-hint-table");
    if (filename == NULL) return;
    char line[1024];
    FILE * f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        /* There's no point in proceeding, since it would really change
         * the behaviour of the program to do so */
        exit(1);
    }
    max_hint_bitsize[0] = 0;
    max_hint_bitsize[1] = 0;
    for(;;) {
        char * x = fgets(line, sizeof(line), f);
        double t;
        unsigned long z;
        /* Tolerate comments and blank lines */
        if (x == NULL) break;
        for( ; *x && isspace(*x) ; x++) ;
        if (*x == '#') continue;
        if (!*x) continue;

        descent_hint h;
        siever_config & sc(h.conf);

        /* start with the global defaults */
        sc = config_base;

        z = strtoul(x, &x, 10);
        ASSERT_ALWAYS(z > 0);
        sc.bitsize = z;
        switch(*x++) {
            case '@' :
                sc.side = strtoul(x, &x, 0);
                ASSERT_ALWAYS(sc.side < 2);
                break;
            default:
                fprintf(stderr, "%s: parse error at %s\n", filename, line);
                exit(1);
        }
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t >= 0);
        h.expected_time = t;
        for( ; *x && isspace(*x) ; x++) ;
        t = strtod(x, &x); ASSERT_ALWAYS(t >= 0);
        h.expected_success = t;
        for( ; *x && isspace(*x) ; x++) ;
        char letter MAYBE_UNUSED = *x;
        for( ; *x && !isdigit(*x) ; x++) ;
        z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
        if (letter == 'I') {
            sc.logI_adjusted = z;
            sc.logA = 2*z-1;
        } else if (letter == 'A') {
            sc.logA = z;
            sc.logI_adjusted = (z+1)/2;
        } else {
            fprintf(stderr, "%s: parse error (want I= or A=) at %s\n", filename, line);
            exit(EXIT_FAILURE);
        }
        
        for(int s = 0 ; s < 2 ; s++) {
            for( ; *x && isspace(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc.sides[s].lim = z;
            for( ; *x && !isdigit(*x) ; x++) ;
            z = strtoul(x, &x, 10); ASSERT_ALWAYS(z > 0);
            sc.sides[s].lpb = z;
            /* recognize this as a double. If it's < 10, we'll consider
             * this means lambda */
            {
                for( ; *x && !isdigit(*x) ; x++) ;
                double t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                if (t < 10) {
                    sc.sides[s].lambda = t;
                    sc.sides[s].mfb = t * sc.sides[s].lpb;
                    /* Then no "lambda" is allowed */
                    continue;
                } else {
                    sc.sides[s].mfb = t;
                }
            }
            if (*x == ',') {
                for( ; *x && !isdigit(*x) ; x++) ;
                t = strtod(x, &x); ASSERT_ALWAYS(t > 0);
                sc.sides[s].lambda = t;
            } else {
                /* this means "automatic" */
                sc.sides[s].lambda = 0;
            }
        }
        for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));

        if (sc.bitsize > max_hint_bitsize[sc.side])
            max_hint_bitsize[sc.side] = sc.bitsize;

        hint_table.push_back(h);
    }
    if (hint_table.empty()) {
        fprintf(stderr, "%s: no data ??\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Allocate the quick lookup tables */

    for(int s = 0 ; s < 2 ; s++) {
        unsigned int n = max_hint_bitsize[s] + 1;
        hint_lookups[s] = (int *)malloc(n * sizeof(int));
        for(unsigned int i = 0 ; i < n ; i++) {
            hint_lookups[s][i] = -1;
        }
    }
    for(unsigned int i = 0 ; i < hint_table.size() ; i++) {
        descent_hint & h(hint_table[i]);
        siever_config & sc = h.conf;
        int d = hint_lookups[sc.side][sc.bitsize];
        if (d >= 0) {
            fprintf(stderr, "Error: two hints found for %d@%d\n",
                    sc.bitsize, sc.side);
            exit(1);
        }
        hint_lookups[sc.side][sc.bitsize] = i;
        /* We must make sure that the default factor base bounds are
         * larger than what we have for all the hint cases, so that it
         * remains reasonable to base our work on the larger factor base
         * (thus doing incomplete sieving).
         */
        /* But now, the .poly file does not contain lim data anymore.
         * This is up to the user to do this check, anyway, because
         * it is not sure that makefb was run with the lim given in the
         * poly file.
         */
        /*
        for(int s = 0 ; s < 2 ; s++) {
            ASSERT_ALWAYS(sc.sides[s]->lim <= cpoly->pols[s]->lim);
        }
        */
    }

    fclose(f);
}/*}}}*/

void las_info::clear_hint_table()
{
    free(hint_lookups[0]);
    free(hint_lookups[1]);
}

#endif /* DLP_DESCENT */

/* {{{ Parse default siever config (fill all possible fields). Return
 * true if the parsed siever config is complete and can be used without
 * per-special-q info. */
bool parse_default_siever_config(siever_config & sc, param_list_ptr pl)
{
    /* The default config is not necessarily a complete bit of
     * information.
     *
     * The field with the bitsize of q is here filled with q0 as found in
     * the command line. Note though that this is only one choice among
     * several possible (q1, or just no default).
     * For the descent, we do not intend to have a default config, thus
     * specifying q0 makes no sense. Likewise, for file-based todo lists,
     * we have no default either, and no siever is configured to provide
     * this ``default'' behaviour (yet, the default bits here are used to
     * pre-fill the config data later on).
     */
    mpz_t q0;
    mpz_init(q0);
    if (param_list_parse_mpz(pl, "q0", q0)) {
        sc.bitsize = mpz_sizeinbase(q0, 2);
    }
    mpz_clear(q0);
    sc.side = 1; // Legacy.
    if (!param_list_parse_int(pl, "sqside", &sc.side)) {
        verbose_output_print(0, 1, "# Warning: sqside not given, "
                "assuming side 1 for backward compatibility.\n");
    }
    param_list_parse_double(pl, "lambda0", &(sc.sides[0].lambda));
    param_list_parse_double(pl, "lambda1", &(sc.sides[1].lambda));
    int complete = 1;
    complete &= param_list_parse_ulong(pl, "lim0", &(sc.sides[0].lim));
    complete &= param_list_parse_ulong(pl, "lim1", &(sc.sides[1].lim));
    if (param_list_lookup_string(pl, "A")) {
        complete &= param_list_parse_int  (pl, "A",    &(sc.logA));
        if (param_list_lookup_string(pl, "I")) {
            fprintf(stderr, "# -A and -I are incompatible\n");
            exit(EXIT_FAILURE);
        }
        // This is useful for a default config
        sc.logI_adjusted = (sc.logA+1)/2;
    } else if (param_list_lookup_string(pl, "I")) {
        int I;
        complete &= param_list_parse_int  (pl, "I", &I);
        sc.logA = 2 * I - 1;
        verbose_output_print(0, 1, "# Interpreting -I %d as meaning -A %d\n", I, sc.logA);
        // This is useful for a default config
        sc.logI_adjusted = I;
    }

    if (sc.sides[0].lim > 2147483647UL || sc.sides[1].lim > 2147483647UL)
    {
        fprintf (stderr, "Error, lim0/lim1 must be < 2^31 (bug #21094)\n");
        exit (EXIT_FAILURE);
    }

    complete &= param_list_parse_int(pl, "lpb0",  &(sc.sides[0].lpb));
    complete &= param_list_parse_int(pl, "mfb0",  &(sc.sides[0].mfb));
    complete &= param_list_parse_int(pl, "lpb1",  &(sc.sides[1].lpb));
    complete &= param_list_parse_int(pl, "mfb1",  &(sc.sides[1].mfb));
    if (!complete) {
        verbose_output_print(0, 1, "# default siever configuration is incomplete ; required parameters are I, lim[01], lpb[01], mfb[01]\n");

    }

    // Sublattices?
    sc.sublat.m = 0; // no sublattices by default.
    param_list_parse_uint(pl, "sublat", &(sc.sublat.m));

    /* Parse optional siever configuration parameters */
    sc.td_thresh = 1024;	/* default value */
    sc.skipped = 1;	/* default value */
    param_list_parse_uint(pl, "tdthresh", &(sc.td_thresh));
    param_list_parse_uint(pl, "skipped", &(sc.skipped));

    sc.unsieve_thresh = 100;
    if (param_list_parse_uint(pl, "unsievethresh", &(sc.unsieve_thresh))) {
        verbose_output_print(0, 1, "# Un-sieving primes > %u\n",
                sc.unsieve_thresh);
    }

    // XXX note that when the sieving range size varies with the
    // special-q, we need to accept that the bucket threshold varies,
    // too.
    //
    // As a consequence, we should make it possible to specify extra
    // parameters in the hint file format (not just I,lim,lpb,mfb).
    sc.bucket_thresh = 1 << ((sc.logA+1)/2);	/* default value */
    sc.bucket_thresh1 = 0;	/* default value */
    sc.bk_multiplier = 0.0;
    /* overrides default only if parameter is given */
    param_list_parse_ulong(pl, "bkthresh", &(sc.bucket_thresh));
    param_list_parse_ulong(pl, "bkthresh1", &(sc.bucket_thresh1));
    param_list_parse_double(pl, "bkmult", &(sc.bk_multiplier));

    const char *powlim_params[2] = {"powlim0", "powlim1"};
    for (int side = 0; side < 2; side++) {
        if (!param_list_parse_ulong(pl, powlim_params[side],
                    &sc.sides[side].powlim)) {
            sc.sides[side].powlim = sc.bucket_thresh - 1;
            verbose_output_print(0, 1,
                    "# Using default value of %lu for -%s\n",
                    sc.sides[side].powlim, powlim_params[side]);
        }
    }

    const char *ncurves_params[2] = {"ncurves0", "ncurves1"};
    for (int side = 0; side < 2; side++)
        if (!param_list_parse_int(pl, ncurves_params[side],
                    &sc.sides[side].ncurves))
            sc.sides[side].ncurves = -1;

    long dupqmin[2] = {0, 0};
    param_list_parse_long_and_long(pl, "dup-qmin", dupqmin, ",");
    sc.sides[0].qmin = dupqmin[0];
    sc.sides[1].qmin = dupqmin[1];

    long dupqmax[2] = {LONG_MAX, LONG_MAX};
    param_list_parse_long_and_long(pl, "dup-qmax", dupqmax, ",");
    sc.sides[0].qmax = dupqmax[0];
    sc.sides[1].qmax = dupqmax[1];

    /* Change qmin = 0 (not initialized) into LONG_MAX */
    for (int side = 0; side < 2; side ++)
      if (sc.sides[side].qmin == 0)
	sc.sides[side].qmin = LONG_MAX;

    /* On the special-q side, if qmin is not given, use
     *   qmin = lim
     * by default.
     */
    if (sc.sides[sc.side].qmin == LONG_MAX)
      sc.sides[sc.side].qmin = sc.sides[sc.side].lim;

    return complete;
}
/* }}} */

las_info::las_info(param_list_ptr pl)/*{{{*/
{
    /* We strive to initialize things in the exact order they're written
     * in the struct */
    // ----- general operational flags {{{
    nb_threads = 1;		/* default value */
    param_list_parse_int(pl, "t", &nb_threads);
    if (nb_threads <= 0) {
	fprintf(stderr,
		"Error, please provide a positive number of threads\n");
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    output = stdout;
    outputname = param_list_lookup_string(pl, "out");
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }
    setvbuf(output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */

    galois = param_list_lookup_string(pl, "galois");
    verbose = param_list_parse_switch(pl, "-v");
    suppress_duplicates = param_list_parse_switch(pl, "-dup");

    verbose_interpret_parameters(pl);
    verbose_output_init(NR_CHANNELS);
    verbose_output_add(0, output, verbose + 1);
    verbose_output_add(1, stderr, 1);
    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, output, 1);
    if (param_list_parse_switch(pl, "-stats-stderr")) {
        /* If we should also print stats to stderr, add stderr to channel 2 */
        verbose_output_add(2, stderr, 1);
    }
#ifdef TRACE_K
    const char *trace_file_name = param_list_lookup_string(pl, "traceout");
    FILE *trace_file = stderr;
    if (trace_file_name != NULL) {
        trace_file = fopen(trace_file_name, "w");
        DIE_ERRNO_DIAG(trace_file == NULL, "fopen", trace_file_name);
    }
    verbose_output_add(TRACE_CHANNEL, trace_file, 1);
#endif
    param_list_print_command_line(output, pl);
    las_display_config_flags();
    /*  Parse polynomial */
    cado_poly_init(cpoly);
    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
	cado_poly_clear(cpoly);
	param_list_clear(pl);
        exit(EXIT_FAILURE);
    }
    if (!cado_poly_read(cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	cado_poly_clear(cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    // sc.skewness = cpoly->skew;
    /* -skew (or -S) may override (or set) the skewness given in the
     * polynomial file */
    param_list_parse_double(pl, "skew", &(cpoly->skew));
    if (cpoly->skew <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	cado_poly_clear(cpoly);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    gmp_randinit_default(rstate);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed))
        gmp_randseed_ui(rstate, seed);
    // }}}

    // ----- default config and adaptive configs {{{
    // perhaps the set of parameters passed is not complete.
    if (parse_default_siever_config(config_base, pl)) {
        default_config_ptr = &config_base;
        sievers.push_back(sieve_info(*this, config_base, pl));
    }
    // }}}

    // ----- stuff roughly related to the descent {{{
#ifdef DLP_DESCENT
    int has_descent_hint = param_list_lookup_string(pl, "descent-hint-table") != NULL;
#else
    const int has_descent_hint = 0;
#endif

    if (!default_config_ptr && !has_descent_hint) {
        fprintf(stderr,
                "Error: no default config set, and no hint table either\n");
        exit(EXIT_FAILURE);
    }

    memset(max_hint_bitsize, 0, sizeof(max_hint_bitsize));
    memset(hint_lookups, 0, sizeof(hint_lookups));
    descent_helper = NULL;
#ifdef  DLP_DESCENT
    init_hint_table(pl);
    // the I value in the default siever config must be at least the max
    // of the I's found in the hint table, because it creates the
    // fb_parts. Let's just gently abort in that case.
    for(unsigned int i = 0 ; i < hint_table.size() ; i++) {
        descent_hint & h(hint_table[i]);
        int logI_adjusted = h.conf.logI_adjusted;
        if (logI_adjusted > config_base.logI_adjusted) {
            fprintf(stderr, "Error: the I value passed in command-line must be at least"
                    " the ones found in the hint file\n");
            exit (EXIT_FAILURE);
        }
    }

    dlog_base = new las_dlog_base();
    dlog_base->lookup_parameters(pl);
    dlog_base->read();
#endif
    // }}}

    // ----- todo list and such {{{
    nq_pushed = 0;
    nq_max = UINT_MAX;
    random_sampling = 0;
    if (param_list_parse_uint(pl, "random-sample", &nq_max)) {
        random_sampling = 1;
    } else if (param_list_parse_uint(pl, "nq", &nq_max)) {
        if (param_list_lookup_string(pl, "q1") || param_list_lookup_string(pl, "rho")) {
            fprintf(stderr, "Error: argument -nq is incompatible with -q1 or -rho\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        todo_list_fd = fopen(filename, "r");
        if (todo_list_fd == NULL) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            cado_poly_clear(cpoly);
            param_list_clear(pl);
            exit(EXIT_FAILURE);
        }
    } else {
        todo_list_fd = NULL;
    }
    // }}}

    // ----- batch mode {{{
    batch = param_list_parse_switch(pl, "-batch");
    batch_print_survivors = param_list_parse_switch(pl, "-batch-print-survivors");
    cofac_list_init (L);
    // }}} 
    
    init_cof_stats(pl);
}/*}}}*/

las_info::~las_info()/*{{{*/
{
    // ----- general operational flags {{{
    if (outputname)
        fclose_maybe_compressed(output, outputname);
    gmp_randclear(rstate);
    verbose_output_clear();
    cado_poly_clear(cpoly);
    // }}}

    // ----- default config and adaptive configs: nothing

    // ----- stuff roughly related to the descent {{{
#ifdef  DLP_DESCENT
    clear_hint_table();
    delete dlog_base;
#endif
    // }}}
 
    // ----- todo list and such {{{
    if (todo_list_fd) {
        fclose(todo_list_fd);
        todo_list_fd = NULL;
    }
    // }}}
 
    // ----- batch mode: very little
    cofac_list_clear (L);

    clear_cof_stats();
}/*}}}*/

// {{{ las_info::{init,clear,print}_cof_stats
void las_info::init_cof_stats(param_list_ptr pl)
{
    const char *statsfilename = param_list_lookup_string (pl, "stats-cofact");
    if (statsfilename != NULL) { /* a file was given */
        if (default_config_ptr == NULL) {
            fprintf(stderr, "Error: option stats-cofact works only "
                    "with a default config\n");
            exit(EXIT_FAILURE);
#ifdef DLP_DESCENT
        } else if (param_list_lookup_string(pl, "descent-hint-table")) {
            verbose_output_print(0, 1, "# Warning: option stats-cofact "
                    "only applies to the default siever config\n");
#endif
        }
        siever_config const & sc(*default_config_ptr);

        cof_stats_file = fopen (statsfilename, "w");
        if (cof_stats_file == NULL)
        {
            fprintf (stderr, "Error, cannot create file %s\n",
                    statsfilename);
            exit (EXIT_FAILURE);
        }
        //allocate cof_call and cof_succ
        int mfb0 = sc.sides[0].mfb;
        int mfb1 = sc.sides[1].mfb;
        cof_call = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        cof_succ = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        for (int i = 0; i <= mfb0; i++) {
            cof_call[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            cof_succ[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            for (int j = 0; j <= mfb1; j++)
                cof_call[i][j] = cof_succ[i][j] = 0;
        }
    } else {
        cof_stats_file = NULL;
        cof_call = NULL;
        cof_succ = NULL;
    }
}
void las_info::print_cof_stats()
{
    if (!cof_stats_file) return;
    ASSERT_ALWAYS(default_config_ptr);
    siever_config const & sc0(*default_config_ptr);
    int mfb0 = sc0.sides[0].mfb;
    int mfb1 = sc0.sides[1].mfb;
    for (int i = 0; i <= mfb0; i++) {
        for (int j = 0; j <= mfb1; j++)
            fprintf (cof_stats_file, "%u %u %u %u\n", i, j,
                    cof_call[i][j],
                    cof_succ[i][j]);
    }
}

void las_info::clear_cof_stats()
{
    if (!cof_stats_file) return;
    ASSERT_ALWAYS(default_config_ptr);
    siever_config const & sc0(*default_config_ptr);
    for (int i = 0; i <= sc0.sides[0].mfb; i++) {
        free (cof_call[i]);
        free (cof_succ[i]);
    }
    free (cof_call);
    free (cof_succ);
    fclose (cof_stats_file);
    cof_stats_file = NULL;
}
//}}}


/* Look for an existing sieve_info in las.sievers with configuration matching
   that in sc; if none exists, create one. */
sieve_info & get_sieve_info_from_config(las_info & las, siever_config const & sc, param_list pl)/*{{{*/
{
    std::list<sieve_info>::iterator psi;
#if 0
#if 0
    psi = find_if(las.sievers.begin(), las.sievers.end(), sc.same_config_q_A_logI());
#else
    psi = find_if(las.sievers.begin(), las.sievers.end(), sc.same_config_q_logI());
#endif
#endif
    psi = find_if(las.sievers.begin(), las.sievers.end(), sc.same_config());
    if (psi != las.sievers.end()) {
        sc.display();
        return *psi;
    }
    las.sievers.push_back(sieve_info(las, sc, pl));
    sieve_info & si(las.sievers.back());
    verbose_output_print(0, 1, "# Creating new sieve configuration for q~2^%d on side %d (logI=%d)\n",
            sc.bitsize, sc.side, si.conf.logI_adjusted);
    sc.display();
    return las.sievers.back();
}/*}}}*/

