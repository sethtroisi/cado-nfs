#include "gmp.h"
#include "utils.h"

typedef struct {
  double * sarr;
  unsigned long B;
  long U0,U1,V0,V1;  
  poly_t f,g,fdg_gdf,gmodp;
  long space;
  int m_max;
  int maxpow;
  mpz_t * ppow;
  double cutoff;
  long l0;
  long igl;
} rootsieve_dictionary_struct;

typedef rootsieve_dictionary_struct rootsieve_dictionary[1];

void rootsieve(poly_t h, poly_t f,poly_t g, unsigned long prime_bound );
