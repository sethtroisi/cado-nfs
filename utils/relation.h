#ifndef RELATION_H_
#define RELATION_H_

#ifndef __cplusplus
#error "This is C++-only"
#endif

#include <stdint.h>
#include <stdio.h>
#include <vector>

#include "cado_poly.h"

struct relation {
    /* Note that for the rational side, we do not compute r !!! */
    struct pr {
        mpz_t p,r;
        int e;
        /* we don't want to use mpz_class, so this is a little bit touchy
         */
        pr() {
            mpz_init(p);
            mpz_init(r);
            e=0;
        }
        pr(mpz_srcptr ap, mpz_srcptr ar, int ae=1) {
            mpz_init_set(p, ap);
            if (ar)
                mpz_init_set(r, ar);
            else
                mpz_init_set_ui(r, 0);
            e = ae;
        }
        pr(unsigned long ap, unsigned long ar, int ae=1) {
            mpz_init_set_ui(p, ap);
            mpz_init_set_ui(r, ar);
            e = ae;
        }
        ~pr() { mpz_clear(p);
            mpz_clear(r);
        }
        pr(pr const& o) { mpz_init_set(p, o.p);
            mpz_init_set(r, o.r);
            e=o.e;
        }
#if __cplusplus >= 201103L
        pr(pr && o) {
            mpz_init(p);
            mpz_init(r);
            mpz_swap(p, o.p);
            mpz_swap(r, o.r);
            e=o.e;
        }
#endif
        pr& operator=(const pr& o) {
            mpz_set(p, o.p);
            mpz_set(r, o.r);
            e=o.e;
            return *this;
        }
    };
    int64_t a;		/* only a is allowed to be negative */
    uint64_t b;
    int rational_side;    /* index of the rational side, if any */
    std::vector<pr> sides[2];

    relation() { a=b=0; }
    operator bool() const { return a || b; }
    relation(int64_t a, uint64_t b, int rational_side = RATIONAL_SIDE)
        : a(a), b(b), rational_side(rational_side)
    {}

    void add(int side, mpz_srcptr p, mpz_srcptr r) {
        sides[side].push_back(pr(p, r));
    }
    void add(int side, unsigned long p, unsigned long r) {
        sides[side].push_back(pr(p, r));
    }

    /* the single-member add() functions recompute r for the algebraic
     * side, based on a,b
     */
    void add(int side, mpz_srcptr p);
    void add(int side, unsigned long p);

    void print(FILE * f, const char * prefix) const;
    int parse(const char * line);


    void fixup_r(bool also_rational = false);
    void compress();
};

#endif	/* RELATION_H_ */
