/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_vbase.h"
#include "abase_u64k1.h"
#include "abase_u64k2.h"
#include "abase_u64k1_t.h"
#include "abase_u64k2_t.h"
void abase_vbase_oo_field_init_byname(abase_vbase_ptr v, const char * s0)
{
    if (strcmp(s0, "u64k1") == 0) {
        abase_u64k1_oo_field_init(v);
    } else if (strcmp(s0, "u64k2") == 0) {
        abase_u64k2_oo_field_init(v);
    } else {
        abort();
    }
}

void abase_vbase_oo_field_init_bygroupsize(abase_vbase_ptr v, int g)
{
    if (g == 64) { abase_vbase_oo_field_init_byname(v, "u64k1"); return; }
    if (g == 128) { abase_vbase_oo_field_init_byname(v, "u64k2"); return; }
    abort();
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
    } else {
        abort();
    }
}


/* vim:set ft=cpp: */
