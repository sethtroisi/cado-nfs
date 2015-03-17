#ifndef RELATION_H_
#define RELATION_H_

#ifndef __cplusplus
#error "This is C++-only"
#endif

#include <stdint.h>
#include <stdio.h>
#include <vector>

#include "cado_poly.h"

struct relation_ab {
    int64_t a;		/* only a is allowed to be negative */
    uint64_t b;
    operator bool() const { return a || b; }
    relation_ab() { a=b=0; }
    relation_ab(int64_t a, uint64_t b) : a(a), b(b) {}
    bool operator<(const relation_ab& o) const {
        return a < o.a || (a == o.a && b < o.b);
    }
};

struct relation : public relation_ab {
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
    int rational_side;    /* index of the rational side, if any */
    std::vector<pr> sides[2];

    relation() {}
    operator bool() const { return (bool) (relation_ab) *this; }
    relation(int64_t a, uint64_t b, int rational_side = RATIONAL_SIDE)
        : relation_ab(a,b), rational_side(rational_side)
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

struct pr_cmp {
    bool operator()(relation::pr const& a, relation::pr const& b) const {
        int c = mpz_cmp(a.p, b.p);
        if (c) { return c < 0; }
        c = mpz_cmp(a.r, b.r);
        return c < 0;
    }
};

#endif	/* RELATION_H_ */
