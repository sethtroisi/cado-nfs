#ifndef __POLYFACTOR_H__
#define __POLYFACTOR_H__

typedef struct {
    int n;
    int alloc;
    fppol_t *factors;
} fppol_fact_struct_t;

typedef fppol_fact_struct_t fppol_fact_t[1];
typedef fppol_fact_struct_t * fppol_fact_ptr;

void fppol_fact_init(fppol_fact_ptr F);
void fppol_fact_clear(fppol_fact_ptr F);
void fppol_fact_push(fppol_fact_ptr F, fppol_t p);
int fppol_fact_pop(fppol_ptr p, fppol_fact_ptr F);
void fppol_fact_out(FILE * out, fppol_fact_ptr F);

void fppol_factor(fppol_fact_ptr factors, fppol_t f);
int fppol_is_irreducible(fppol_srcptr f);




#endif   /* __POLYFACTOR_H__ */
