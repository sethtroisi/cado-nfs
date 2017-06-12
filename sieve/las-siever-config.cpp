#include "cado.h"
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include "las-forwardtypes.hpp"
#include "las-siever-config.hpp"
#include "las-types.hpp"
#include "verbose.h"

/* siever_config stuff */

void siever_config::display() const /*{{{*/
{
    if (bitsize == 0) return;

    verbose_output_print(0, 1, "# Sieving parameters for q~2^%d on side %d\n",
            bitsize, side);
    /* Strive to keep these output lines untouched */
    verbose_output_print(0, 1,
	    "# Sieving parameters: lim0=%lu lim1=%lu lpb0=%d lpb1=%d\n",
	    sides[0].lim, sides[1].lim,
            sides[0].lpb, sides[1].lpb);
    verbose_output_print(0, 1,
	    "#                     mfb0=%d mfb1=%d\n",
	    sides[0].mfb, sides[1].mfb);
    if (sides[0].lambda != 0 || sides[1].lambda != 0) {
        verbose_output_print(0, 1,
                "#                     lambda0=%1.1f lambda1=%1.1f\n",
            sides[0].lambda, sides[1].lambda);
    }
}/*}}}*/

siever_config las_info::get_config_for_q(las_todo_entry const & doing) const /*{{{*/
{
    // arrange so that we don't have the same header line as the one
    // which prints the q-lattice basis
    verbose_output_vfprint(0, 1, gmp_vfprintf,
                         "#\n"
                         "# "
                         "Now sieving side-%d q=%Zd; rho=%Zd\n",
                         doing.side,
                         (mpz_srcptr) doing.p,
                         (mpz_srcptr) doing.r);
    siever_config config = config_base;
    config.bitsize = mpz_sizeinbase(doing.p, 2);
    config.side = doing.side;

    /* Do we have a hint table with specifically tuned parameters,
     * well suited to this problem size ? */
    for(unsigned int i = 0 ; i < hint_table.size() ; i++) {
        siever_config const & sc(hint_table[i].conf);
        if (!sc.same_config_q()(config)) continue;
        verbose_output_print(0, 1, "# Using parameters from hint list for q~2^%d on side %d [%d@%d]\n", sc.bitsize, sc.side, sc.bitsize, sc.side);
        config = sc;
    }

    if (doing.iteration) {
        verbose_output_print(0, 1, "#\n# NOTE:"
                " we are re-playing this special-q because of"
                " %d previous failed attempt(s)\n", doing.iteration);
        /* update sieving parameters here */
        double ratio = double(config.sides[0].mfb) /
            double(config.sides[0].lpb);
        config.sides[0].lpb += doing.iteration;
        config.sides[0].mfb = ratio*config.sides[0].lpb;
        ratio = double(config.sides[1].mfb) /
            double(config.sides[1].lpb);
        config.sides[1].lpb += doing.iteration;
        config.sides[1].mfb = ratio*config.sides[1].lpb;
        verbose_output_print(0, 1,
                "# NOTE: current values of lpb/mfb: %d,%d %d,%d\n#\n", 
                config.sides[0].lpb,
                config.sides[0].mfb,
                config.sides[1].lpb,
                config.sides[1].mfb);
    }

    return config;
}/*}}}*/

