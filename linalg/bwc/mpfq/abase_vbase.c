/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include <stdarg.h>
#include "abase_vbase.h"
/* hello world ! */
#include "abase_u64k1.h"
#include "abase_u64k1_t.h"
#include "abase_u64k2.h"
#include "abase_u64k2_t.h"
#include "abase_p16.h"
#include "abase_p16_t.h"
#include "abase_p_4.h"
#include "abase_p_4_t.h"
void abase_vbase_oo_field_init_byfeatures(abase_vbase_ptr v, ...)
{
        va_list ap;
        va_start(ap, v);
        mpz_t p;
        mpz_init_set_ui(p, 2);
        int groupsize = 1;
        for(int a ; (a = va_arg(ap, int)) != 0 ; ) {
            if (a == MPFQ_PRIME_MPZ) {
                mpz_set(p, va_arg(ap, mpz_srcptr));
            } else if (a == MPFQ_GROUPSIZE) {
                groupsize = va_arg(ap, int);
            } else {
                /* We do not support MPFQ_PRIME_MPN. Only MPFQ_PRIME_MPZ*/
                fprintf(stderr, "Feature code %d unsupported\n", a);
                exit(1);
            }
        }
        va_end(ap);
        if (mpz_cmp_ui(p, 2) == 0) {
            if (groupsize == 64) {
                abase_u64k1_oo_field_init(v);
            } else if (groupsize == 128) {
                /* TODO: SSE-2 ! */
                abase_u64k2_oo_field_init(v);
            } else {
                fprintf(stderr, "Unsupported group size %d\n", groupsize);
                exit(1);
            }
        } else if (groupsize > 1) {
            /* All prime fields require group size == 1 */
            fprintf(stderr, "Unsupported group size %d\n", groupsize);
            exit(1);
        } else if (mpz_cmp_ui(p, 1 << 16) < 0) {
            abase_p16_oo_field_init(v);
        } else if (mpz_size(p) == 4) { abase_p_4_oo_field_init(v);
        } else {
            fprintf(stderr, "Unsupported prime (%zu limbs)\n", mpz_size(p));
            exit(1);
        }
        v->field_specify(v, MPFQ_PRIME_MPZ, p);
        v->field_specify(v, MPFQ_GROUPSIZE, &groupsize);
        mpz_clear(p);
}

void abase_vbase_oo_init_templates(abase_vbase_tmpl_ptr w, abase_vbase_ptr v0, abase_vbase_ptr v1)
{
    const char * s0 = v0->oo_impl_name(v0);
    const char * s1 = v1->oo_impl_name(v1);
    if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "u64k1") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k1_u64k1_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_u64k1_u64k1_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_u64k1_u64k1_transpose;
    } else if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "u64k2") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k1_u64k2_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_u64k1_u64k2_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_u64k1_u64k2_transpose;
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "u64k1") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_u64k1_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_u64k2_u64k1_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_u64k2_u64k1_transpose;
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "u64k2") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_u64k2_u64k2_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_u64k2_u64k2_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_u64k2_u64k2_transpose;
    } else if (strcmp(s0, "p16") == 0 && strcmp(s1, "p16") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_p16_p16_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_p16_p16_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_p16_p16_transpose;
    } else if (strcmp(s0, "p_4") == 0 && strcmp(s1, "p_4") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_p_4_p_4_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_p_4_p_4_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_p_4_p_4_transpose;
    } else {
        abort();
    }
}


/* vim:set ft=cpp: */
