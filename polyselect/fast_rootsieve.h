typedef struct {
  double * sarr;
  int B;
  double global_offset;
  mpz_poly * f,g,fdg_gdf,gmodp;
  long umax,vmax,space;
  int maxpow;
  long * ppow;
  double cutoff;
  long l0;
  long igl;
} rootsieve_dictionary;

// Can we do this with mod_ul?, directly?
long modular_inverse( long r, long N );

// Rootsieve - Specific
long rotation_inner(rootsieve_dictionary rdict, long p,mpz_poly * ff,mpz_poly * gg,mpz_poly * fdg_gdf,long u0,long v0,long l0,int ld,int m,long twist_v1,double scale,mpz_poly dphi);
int fill_alpha(rootsieve_dictionary rdict, long p);
void extend_ppow_list(rootsieve_dictionary rdict,int l, long p);

// Polynomials
void poly_alloc_identity(poly_t f); 
void compose_psi(poly_t f, poly_t g,unsigned long l,rootsieve_dictionary rdict);
void compose_reduce(poly_t f, poly_t g,unsigned long l,rootsieve_dictionary rdict);
