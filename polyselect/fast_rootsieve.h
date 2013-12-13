#include "gmp.h"
#include "utils.h"

typedef struct {
  // Stocks the lattice of origin (u0,v0) with basis (us,skew),(0,vs), associated to prime p.
  unsigned long p; 
  long u0,v0;
  long us;
  long skew;
  long vs;
  double contribution;
} rootsieve_lattice;

//typedef rootsieve_lattice_struct rootsieve_lattice[1];

typedef struct {
  rootsieve_lattice * lattices;
  long length;
} lattice_list;

typedef struct {
  rootsieve_lattice * lattices;
  long lattice_alloc;
  long L;
  long U0,U1,V0,V1;  
  mpz_poly_t f,g,fdg_gdf,gmodp;
  int m_max;
  int maxpow;
  mpz_t * ppow;
  double cutoff;
  long l0;
  long igl;
} rootsieve_dictionary_struct;

typedef rootsieve_dictionary_struct rootsieve_dictionary[1];

lattice_list rootsieve(mpz_t * f, int degf, mpz_t m, mpz_t b, long U0, long U1, long V0, long V1, unsigned long prime_bound);
void print_lattices(lattice_list l);
