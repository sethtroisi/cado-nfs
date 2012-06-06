/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_vbase.h"
#include "abase_p16.h"
#include "abase_p16_t.h"
void abase_vbase_oo_field_init_byname(abase_vbase_ptr v, const char * s0)
{
    if (strcmp(s0, "p16") == 0) {
        abase_p16_oo_field_init(v);
    } else {
        abort();
    }
}

void abase_vbase_oo_field_init_bygroupsize(abase_vbase_ptr v, int g MAYBE_UNUSED)
{
        assert(g == 1);
        abase_vbase_oo_field_init_byname(v, "p16");
        return;
}

void abase_vbase_oo_init_templates(abase_vbase_tmpl_ptr w, abase_vbase_ptr v0, abase_vbase_ptr v1)
{
    const char * s0 = v0->oo_impl_name(v0);
    const char * s1 = v1->oo_impl_name(v1);
    if (strcmp(s0, "p16") == 0 && strcmp(s1, "p16") == 0) {
        w->dotprod = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int)) abase_p16_p16_dotprod;
        w->addmul_tiny = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int)) abase_p16_p16_addmul_tiny;
        w->transpose = (void (*) (abase_vbase_ptr, abase_vbase_ptr, void *, const void *)) abase_p16_p16_transpose;
    } else {
        abort();
    }
}


/* vim:set ft=cpp: */
