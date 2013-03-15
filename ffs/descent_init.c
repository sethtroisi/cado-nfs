#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "timing.h"
#include "params.h"
#include "smoothness.h"
#include "polyfactor.h"


/* Half extended Euclidian algorithm between p and m.
   We want to stop it when u*m + v*p = r with deg(v) ~ deg(r) ~ deg(m)/2 
 /!\ Assume that m is monic and that p is reduced modulo m. */
int fppol_half_extended_Euclidian(fppol_ptr r, fppol_ptr v, fppol_srcptr p, fppol_srcptr m)
{
  /* We will see later if these requirements are a burden to    */
  /* the caller or not.                                         */
  ASSERT_ALWAYS(fppol_is_monic(m));
  ASSERT_ALWAYS(fppol_deg(p) < fppol_deg(m));
  fppol_t r0, r1, v0, v1, tmp;
  fppol_inits(r0, r1, v0, v1, tmp, NULL);

  /* Euclid-reduce (r0, r1) and maintain the following:         */
  /*   u0*m + v0*p = r0   (u0, u1 not computed!)                */
  /*   u1*m + v1*p = r1                                         */
  fppol_set(r0, m);
  fppol_set(r1, p);
  fppol_set_zero(v0);
  fppol_set_one (v1);
  int halfdm = fppol_deg(m) / 2;
  int d0 = fppol_deg(r0);
  int d1 = fppol_deg(r1);
  IF(CMP(FP_SIZE, 2), EQ, , 
     fp_t lc0; fp_t lc1; fp_t lc;
     fppol_get_coeff(lc0, r0, d0);
     fppol_get_coeff(lc1, r1, d1);
  )
  for (int d = d0 - d1; ; ) {
    fppol_mul_ti(tmp, r1, d);
    IF(CMP(FP_SIZE, 2), EQ, ,
       fp_div(lc, lc0, lc1);
       fppol_smul(tmp, tmp, lc);
    )
    fppol_sub(r0, r0, tmp);

    fppol_mul_ti(tmp, v1, d);
    IF(CMP(FP_SIZE, 2), EQ, ,
       fppol_smul(tmp, tmp, lc);
    )
    fppol_sub(v0, v0, tmp);

    d0 = fppol_deg(r0);
    IF(CMP(FP_SIZE, 2), EQ, ,
       if (d0 >= 0) fppol_get_coeff(lc0, r0, d0);
    )

    if (d0 <= halfdm) {
      IF(CMP(FP_SIZE, 2), EQ,
	 fppol_set (v, v0);
	 fppol_set (r, r0);, 
	 fppol_sdiv(v, v0, lc0);
	 fppol_sdiv(r, r0, lc0);
      )
      fppol_clears(r1, v1, tmp, NULL);
      return d0;
    }

    d = d0 - d1;
    if (d <= 0) {
      fppol_swap(r0, r1);
      fppol_swap(v0, v1);
      IF(CMP(FP_SIZE, 2), EQ, , fp_swap(lc0, lc1);)
      int dt = d0; d0 = d1; d1 = dt;
      d = -d;
    }
  }
  ASSERT_ALWAYS(0); /* Should never go there. */
}


void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  lpb *           large prime bound \n");
    fprintf(stderr, "  target *        element whose discrete log has to be computed \n");
    fprintf(stderr, "  phi *           finite field defining polynomial \n");
    /* fprintf(stderr, "  beta *         generator of the subgroup target belongs to \n"); */
    fprintf(stderr, "  gf              indicates the base field for sanity check\n");
 
    /* We consider we can assume that t is a "valid" generator even if it is not,
     hence beta is commented for the moment */

    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(1);
}

int main(int argc, char **argv)
{
    int lpb = 0; 
    fppol_t target, phi; /*,beta */ 
    char *argv0 = argv[0];
    int gf = 0;

    param_list pl;
    param_list_init(pl);
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Could also be a parameter file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL);
    }
    // Update parameter list at least once to register argc/argv pointers.
    param_list_update_cmdline(pl, &argc, &argv);
    param_list_print_command_line(stdout, pl);

    param_list_parse_int(pl, "gf", &gf);
    if (gf) {
        if (gf != FP_SIZE) {
            fprintf(stderr, "Error: base field mismatch.\n");
            fprintf(stderr, "  The binary is compiled for GF(%d)\n", FP_SIZE);
            fprintf(stderr, "  The parameters are for GF(%d)\n", gf);
            exit(EXIT_FAILURE);
        }
    }
    // read various bounds
    param_list_parse_int(pl, "lpb", &lpb);
    if (lpb == 0) usage(argv0, "lpb");

    // read function field polynomials
    {
        const char * polstr;
        fppol_init(target);
        fppol_init(phi);
        /* fppol_init(beta); */
        polstr = param_list_lookup_string(pl, "target");
        if (polstr == NULL) usage(argv0, "target");
        fppol_set_str(target, polstr);
        polstr = param_list_lookup_string(pl, "phi");
        if (polstr == NULL) usage(argv0, "phi");
        fppol_set_str(phi, polstr);
        /* polstr = param_list_lookup_string(pl, "beta");
        if (polstr == NULL) usage(argv0, "beta");
        fppol_set_str(beta, polstr); */
    }

    // Reduce target modulo phi
    fppol_rem(target, target, phi);

    /* Randomizing target : target = target * t^j
       for j = 1, 2, 3... 
     */
    fppol_t r, v;
    fppol_inits(r, v, NULL);
    unsigned int j;
    for (j = 0; ; ++j) {
      /* Applying an extended Euclidian algoritm to target and phi
       to stop when u.phi + v. target = r whit deg(v) ~ deg(r) ~ n/2 */
      fppol_half_extended_Euclidian(r, v, target, phi);
	
      /* Apply a smoothness test on v and r to check if they are lpb-smooth */
      if (fppol_is_smooth(r, lpb) && fppol_is_smooth(v, lpb))  
         break;
      
      /* multiply target by t for next try */
      fppol_multmod(target, target, phi);
    }

    /* If v and r are lpb-smooth, we get log target + j = log r - log v 
       which involves only polynomials of degree less than lpb:
       the initialization of the special-q descent is over */

    /* We want to output j, r and v */
    fprintf(stdout, "j = %u\n", j);
    fprintf(stdout, "r = ");
    fppol_out(stdout, r);
    fprintf(stdout, "\n");
	 fprintf(stdout, "v = ");
    fppol_out(stdout, v);
    fprintf(stdout, "\n");
    param_list_clear(pl);

    // TODO: clear all fppol_t at the end...
    
    return EXIT_SUCCESS;
}
