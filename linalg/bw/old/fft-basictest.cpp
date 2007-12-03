#include <cstdlib>
#include <gmp.h>
#include <list>
#include "ops_poly_wrap.hpp"
#include "bw_scalar.h"
#include "field_usage.h"

using namespace std;

typedef ops_poly_wrap ft_order_t;

void mb_transform(ft_order_t const& s,
		ft_order_t::mat_mb p, bw_mbpoly src, int d)
{
    for (int i = 0; i < m_param; i++) { 
        for (int j = 0; j < bigdim; j++) {
            mp_limb_t * src0;
            mp_limb_t * src1;
            src0 = mbmat_scal(mbpoly_coeff(src, 0), i, j);
            src1 = mbmat_scal(mbpoly_coeff(src, 1), i, j);
            s.transform(s.mat_mb_get(p, i, j), src0, src1-src0, d);
    
        }
    }
}

void bb_transform(ft_order_t const& s,
		ft_order_t::mat_bb p, bw_bbpoly src, int d)
{
    for (int i = 0; i < bigdim; i++) { 
        for (int j = 0; j < bigdim; j++) {
            mp_limb_t * src0;
            mp_limb_t * src1;
            src0 = bbmat_scal(bbpoly_coeff(src, 0), i, j);
            src1 = bbmat_scal(bbpoly_coeff(src, 1), i, j);
            s.transform(s.mat_bb_get(p, i, j), src0, src1-src0, d);
    
        }
    }
}

void mb_itransform(ft_order_t const& s,
		bw_mbpoly dest, ft_order_t::mat_mb p, int d)
{
    for (int i = 0; i < m_param; i++) { 
        for (int j = 0; j < bigdim; j++) {
            mp_limb_t * dst0;
            mp_limb_t * dst1;
            dst0 = mbmat_scal(mbpoly_coeff(dest, 0), i, j);
            dst1 = mbmat_scal(mbpoly_coeff(dest, 1), i, j);
            s.itransform(dst0, dst1-dst0, d, s.mat_mb_get(p, i, j));
        }
    }
}

void convolution(ft_order_t const& s,
		ft_order_t::mat_mb c,
		ft_order_t::mat_mb a, ft_order_t::mat_bb b)
{
    for (int i = 0; i < m_param; i++) {
        for (int j = 0; j < bigdim; j++) {
            s.zero(s.mat_mb_get(c, i, j));
            for (int l = 0; l < bigdim; l++) {
                s.convolution(
                        s.mat_mb_get(c, i, j),
                        s.mat_mb_get(a, i, l),
                        s.mat_bb_get(b, l, j));
            }
        }
    }
}


int do_one(const char * prime, int dm, int dn, int db)
{
	std::list<char *> nothing;

	printf("%dx%d * %dx%d matrices mod %s\n", dm, dn, dn, db, prime);
	fflush(stdout);

	bw_read_modulus_info(prime,1);

#ifdef HARDCODE_PARAMS
#error "ouch"
#endif

	computed_m_param = dm;
	computed_n_param = dn;
	computed_bigdim = db;

	ft_order_t::init(10000, nothing);

	ft_order_t s;

	bw_mbpoly a;
	bw_bbpoly b;
	bw_mbpoly c;

	ft_order_t::mat_mb ta, tc;
	ft_order_t::mat_bb tb;
	for(int d = 1 ; d < 40 ; d++) {
		// set takes a number of coefficients
		s.set(d*2+1);
		mbpoly_alloc(a, d);
		bbpoly_alloc(b, d);
		mbpoly_alloc(c, 2*d);

		mbpoly_zero(a, d);
		bbpoly_zero(b, d);
		mbpoly_zero(c, 2*d);

		s.mat_mb_alloc(ta);
		s.mat_bb_alloc(tb);
		s.mat_mb_alloc(tc);
		for(int i = 0 ; i < m_param ; i++) {
			for(int j = 0 ; j < bigdim ; j++) {
				s.zero(s.mat_mb_get(ta,i,j));
				s.zero(s.mat_mb_get(tc,i,j));
			}
		}
		for(int i = 0 ; i < bigdim ; i++) {
			for(int j = 0 ; j < bigdim ; j++) {
				s.zero(s.mat_bb_get(tb,i,j));
			}
		}

		for(int k = 0 ; k <= d ; k++) {
			for(int i = 0 ; i < m_param ; i++) {
				for(int j = 0 ; j < bigdim ; j++) {
					k_set_int(mbmat_scal(mbpoly_coeff(a,k),i,j),rand());
				}
			}
			for(int i = 0 ; i < bigdim ; i++) {
				for(int j = 0 ; j < bigdim ; j++) {
					k_set_int(bbmat_scal(bbpoly_coeff(b,k),i,j),rand());
				}
			}
		}

		mb_transform(s, ta, a, d);
		bb_transform(s, tb, b, d);
		convolution(s, tc, ta, tb);
		mb_itransform(s, c, tc, 2*d);

		s.mat_mb_free(ta);
		s.mat_bb_free(tb);
		s.mat_mb_free(tc);

		mbpoly_free(a);
		bbpoly_free(b);
		mbpoly_free(c);
	}

	ft_order_t::cleanup();

	return 0;
}

int main()
{
	for(int m = 1 ; m < 4 ; m++) {
		for(int n = 1 ; n < 4 ; n++) {
			for(int b = 1 ; b < 4 ; b++) {
				do_one("65537", m, n, b);
			}
		}
	}
	return 0;
}

