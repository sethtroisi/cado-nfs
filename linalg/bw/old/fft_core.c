#include <gmp.h>
#include <stdlib.h>
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
#include "fft_core.h"
#include "manu.h"

mp_limb_t *l_roots;
int *phi_tab;
int tabulated_order;
fft_capable_function direct_fft, inverse_fft;

/**********************************************************************/
/* fft engine. */

/* compute s1 <-- s1 + omega * s2
 *         s2 <-- s1 - omega * s2
 */
void butterfly(mp_limb_t * s1, mp_limb_t * s2, mp_limb_t * omega)
{
    mp_limb_t *tmp;

    tmp = FAST_ALLOC(l_size * sizeof(mp_limb_t));
    l_mul(tmp, s2, omega);
    l_sub(s2, s1, tmp);
    l_add(s1, s1, tmp);
    FAST_FREE(tmp);
}

/* compute s1 <--  s1 + s2
 *         s2 <-- (s1 - s2) * omega
 */
void butterfly_inverse(mp_limb_t * s1, mp_limb_t * s2, mp_limb_t * omega)
{
    mp_limb_t *tmp;

    tmp = FAST_ALLOC(l_size * sizeof(mp_limb_t));
    l_sub(tmp, s1, s2);
    l_add(s1, s1, s2);
    l_mul(s2, tmp, omega);
    FAST_FREE(tmp);
}

#if 0
/*
 * compute s1 <-- s1 + s2
 *         s2 <-- s1 - s2
 */
void basic_butterfly(mp_limb_t * s1, mp_limb_t * s2)
{
    mp_limb_t *tmp;

    tmp = FAST_ALLOC(l_size * sizeof(mp_limb_t));
    l_add(tmp, s1, s2);
    l_sub(s2, s1, s2);
    l_set(s1, tmp);
    FAST_FREE(tmp);
}
#endif


static void compute_phi(int d)
{
    int i, j, l;
    phi_tab[0] = 0;
    for (j = 1; j < (1 << d); phi_tab[j++] = 123124);	/* on fout la merde... */
    for (l = 1; l <= d; l++) {
	for (i = 0; i < (1 << (l - 1)); i++) {
	    phi_tab[(i << (d - l + 1)) + (1 << (d - l))] =
		phi_tab[i << (d - l + 1)] + (1 << (l - 1));
	}
    }
}

static unsigned int t_func(unsigned int n)
{
    unsigned int k;
    if (n == 0)
	return 0;
    for (k = 1; k <= n; k <<= 1);
    return n ^ ((k >> 1) - 1);
}

/* this is the core of the FFT. Is works for any l with the good
 * properties */

void direct_fft_recursive(mp_limb_t * buf, mp_size_t gap,
			  mp_limb_t * roots, unsigned int d)
{
    unsigned int i;
    mp_limb_t *p, *q, *r;
    if (d > 1) {
	direct_fft_recursive(buf, gap << 1, roots, d - 1);
	direct_fft_recursive(buf + gap, gap << 1, roots, d - 1);
    }
    p = buf;
    q = buf + gap;
    r = roots;
    for (i = 0; i < (1 << (d - 1)); i++) {
	butterfly(p, q, r);
	p += gap << 1;
	q += gap << 1;
	r += l_size << 1;
    }
}

void inverse_fft_recursive(mp_limb_t * buf, mp_size_t gap,
			   mp_limb_t * roots, unsigned int d)
{
    unsigned int i;
    mp_limb_t *p, *q;
    p = buf;
    q = buf + gap;
    for (i = 0; i < (1 << (d - 1)); i++) {
	butterfly_inverse(p, q, roots + t_func(i << 1) * l_size);
	p += gap << 1;
	q += gap << 1;
    }
    if (d > 1) {
	inverse_fft_recursive(buf, gap << 1, roots, d - 1);
	inverse_fft_recursive(buf + gap, gap << 1, roots, d - 1);
    }
}

void direct_fft_iterative(mp_limb_t * buf, mp_size_t gap,
			  mp_limb_t * roots, unsigned int d)
{
    unsigned int l, offset, k;
    for (l = 1; l <= d; l++) {
	for (offset = 0; offset < (1 << (d - l)); offset++) {
	    for (k = 0; k < (1 << (l - 1)); k++) {
		butterfly(buf + (offset + (k << (d - l + 1))) * gap,
			  buf + (offset + (((k << 1) + 1) << (d - l)) * gap),
			  roots + (k << 1) * l_size);
	    }
	}
    }
}

void inverse_fft_iterative(mp_limb_t * buf, mp_size_t gap,
			   mp_limb_t * roots, unsigned int d)
{
    unsigned int l, offset, k;
    for (l = d; l >= 1; l--) {
	for (offset = 0; offset < (1 << (d - l)); offset++) {
	    for (k = 0; k < (1 << (l - 1)); k++) {
		butterfly_inverse(buf + (offset + (k << (d - l + 1))) * gap,
				  buf + (offset +
					 (((k << 1) + 1) << (d - l)) * gap),
				  roots + t_func(k << 1) * l_size);
	    }

	}

    }
}

void prepare_fft_engine(unsigned int max_order, int enable_cplx)
{
    mp_limb_t *r;
    int i;
    mp_limb_t *ptr, *nptr;

    tabulated_order = max_order;
    phi_tab = malloc((1 << max_order) * sizeof(int));
    compute_phi(max_order);

    prepare_fields_for_fft(max_order, enable_cplx);
    r = fetch_primitive_root(max_order);

    l_roots = malloc((1 << max_order) * l_size * sizeof(mp_limb_t));
    ptr = l_roots;		/* phi_tab[0]=0. Always. */
    l_set_one(ptr);
    for (i = 1; i < (1 << max_order); i++) {
	nptr = l_roots + phi_tab[i] * l_size;
	l_mul(nptr, ptr, r);
	ptr = nptr;
    }
    free(r);

    /* Just a matter of taste, since the difference in terms of
     * efficiency is hardly noticeable */

    direct_fft = &direct_fft_recursive;
    inverse_fft = &inverse_fft_recursive;
}

void cleanup_fft_engine(void)
{
    free(phi_tab);
    free(l_roots);
}

