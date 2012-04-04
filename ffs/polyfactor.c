#include "fppol.h"
#include "types.h"
#include "polyfactor.h"
#include "fppol_facttools.h"



// For lack of a better solution, use deterministic random: take
// polynomials in lex order (maybe random enough for an edf).
static void fppol_set_random(fppol_ptr f, int deg) {
#ifdef USE_F2
    if (deg == 0)
        return;
    ij_t r;
    r[0] = rand();
    int bits = MIN(deg+1, 16);
    r[0] &= (1u << bits)-1;
    fppol_set_ij(f, r);
#else
    static int first_time = 1;
    static ij_t state;

    if (deg == 0) {
        first_time = 1;
        return;
    }

    if (first_time) {
        first_time = 0;
        ij_set_ti(state, 1);
    } else {
        if (!ij_set_next(state, state, MIN(deg+1, 16)))
            ij_set_ti(state, 1);
    }
    fppol_set_ij(f, state);
#endif
}


void fppol_fact_init(fppol_fact_ptr F)
{
    F->alloc = 4;
    F->factors = (fppol_t *) malloc(F->alloc*sizeof(fppol_t));
    for (int i = 0; i < F->alloc; ++i)
        fppol_init(F->factors[i]);
    F->n = 0;
}

void fppol_fact_clear(fppol_fact_ptr F)
{
    for (int i = 0; i < F->alloc; ++i)
        fppol_clear(F->factors[i]);
    free(F->factors);
}

void fppol_fact_push(fppol_fact_ptr F, fppol_t p)
{
    if (F->n == F->alloc) {
        int margin = 4;
        F->factors = (fppol_t *) realloc(F->factors,
                (F->alloc+margin)*sizeof(fppol_t));
        for (int i = 0; i < margin; ++i)
            fppol_init(F->factors[F->alloc + i]);
        F->alloc += margin;
    }
    fppol_set(F->factors[F->n], p);
    F->n++;
}


int fppol_fact_pop(fppol_ptr p, fppol_fact_ptr F)
{
    if (F->n == 0)
        return 0;
    fppol_set(p, F->factors[F->n-1]);
    F->n--;
    return 1;
}

void fppol_fact_out(FILE *out, fppol_fact_ptr F)
{
    for (int i = 0; i < F->n; ++i) {
        fppol_out(out, F->factors[i]);
        if (i != F->n - 1)
            fprintf(out, ",");
    }
}




static void fppol_edf(fppol_fact_ptr factors, fppol_t f, int d)
{
    fppol_t a, tra;
    if (fppol_deg(f) == d) {
        fppol_fact_push(factors, f);
        return;
    }
    ASSERT(fppol_deg(f) % d == 0);

    fppol_init(a);
    fppol_init(tra);

    do {
        fppol_set_random(a, fppol_deg(f)-1); // like in textbooks !
        fppol_set_zero(tra);
        for (int i = 0; i < d-1; ++i) {
            fppol_add(tra, tra, a);
            fppol_qpow(a, a);
            fppol_rem(a, a, f);
        }
        fppol_add(tra, tra, a);

        fppol_gcd(tra, tra, f);
    } while (fppol_deg(tra) == 0 || fppol_deg(tra) == fppol_deg(f));

    ASSERT(fppol_deg(tra) % d == 0);

    if (fppol_deg(tra) == d) {
        fppol_fact_push(factors, tra);
    } else {
        fppol_edf(factors, tra, d);
    }

    fppol_div(tra, f, tra);

    if (fppol_deg(tra) == d) {
        fppol_fact_push(factors, tra);
    } else {
        fppol_edf(factors, tra, d);
    }
    fppol_clear(a);
    fppol_clear(tra);
}


// Trial divide F with the factors given in fact.
// They end-up (with multiplicity) in factors.
static void purge_factors(fppol_ptr F, fppol_fact_ptr factors,
        fppol_fact_ptr fact)
{
    fppol_t quo, rem;
    fppol_init(quo);
    fppol_init(rem);
    for (int i = 0; i < fact->n; ++i) {
        do {
            fppol_divrem(quo, rem, F, fact->factors[i]);
            if (fppol_is_zero(rem)) {
                fppol_set(F, quo);
                fppol_fact_push(factors, fact->factors[i]);
            }
        } while (fppol_is_zero(rem));
    }
    fppol_clear(quo);
    fppol_clear(rem);
}



// Factorize into *monic* irreducible polynomials.
// If F is not monic, it is transformed to monic before factorization.
void fppol_factor(fppol_fact_ptr factors, fppol_t f)
{
    int i;

    fppol_t F;
    fppol_t tqi;
    fppol_t t;
    fppol_t g;

    fppol_init(F);
    fppol_init(t);
    fppol_init(tqi);
    fppol_init(g);

    fppol_fact_t fact;
    fppol_fact_init(fact);
    factors->n = 0;  // start with a fresh data structure (init done by caller)

    // Work with a (monic) copy of f.
    {
        fp_t lc;
        fppol_get_coeff(lc, f, fppol_deg(f));
        fppol_sdiv(F, f, lc);
    }

    i = 0;
    fppol_set_ti(t, 1);
    fppol_set(tqi, t);

    while (fppol_deg(F) > 0) {
        i++;
        fppol_qpow(tqi, tqi);
        fppol_rem(tqi, tqi, F); 
        fppol_sub(g, tqi, t);
        fppol_gcd(g, g, F);
        if (fppol_deg(g) > 0) {
            fppol_set_random(NULL, 0); // reinit the pseudo-random: TODO: Beurk!
            fppol_edf(fact, g, i);
            purge_factors(F, factors, fact);
            fact->n = 0;
        }
    }
    
    fppol_clear(F);
    fppol_clear(t);
    fppol_clear(tqi);
    fppol_clear(g);

    fppol_fact_clear(fact);
}
