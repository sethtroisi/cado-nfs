#include "gmp.h"
#include "utils.h"

typedef struct {
  double * sarr;
  int B;
  double global_offset;
  poly_t f,g,fdg_gdf,gmodp;
  long umax,vmax,space;
  int maxpow;
  mpz_t * ppow;
  double cutoff;
  long l0;
  long igl;
} rootsieve_dictionary_struct;

typedef rootsieve_dictionary_struct rootsieve_dictionary[1];
// We can eliminate the following by using mod_ul residues directly
unsigned long modular_inverse( unsigned long r, unsigned long N );
unsigned long modular_product( unsigned long a, unsigned long b, unsigned long N ); 
// Rootsieve - Specific
void rootsieve(poly_t h, poly_t f,poly_t g, long U, long V, unsigned long prime_bound );
long rotation_inner(rootsieve_dictionary rdict, unsigned long p,poly_t ff,poly_t gg,poly_t fdg_gdf,long u0,long v0,long l0,int ld,int m,long twist_v1,double scale,poly_t dphi);
long fill_alpha(rootsieve_dictionary rdict, unsigned long p);
long light (double * alpha, long u0, long v0, long us, long vs, long skew, long um, long vm, double contribution);
void extend_ppow_list(rootsieve_dictionary rdict,int newlength, unsigned long p); 

// Polynomials
void poly_alloc_identity(poly_t f); 
void compose_psi(poly_t f, const poly_t g,unsigned long l,rootsieve_dictionary rdict);
void compose_reduce(poly_t f, const poly_t g,unsigned long l,rootsieve_dictionary rdict, mpz_t m);
