typedef struct {
  int deg;
  mpz_t * f;
} mpz_poly;

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

void mpz_poly_init(mpz_poly * p, int deg);

void mpz_poly_clear(mpz_poly * p);

void mpz_poly_trim(mpz_poly * p);

mpz_poly mpz_poly_derivative(mpz_poly * f);

mpz_poly mpz_poly_coeff_reduction(mpz_poly * f, unsigned long m);

mpz_poly mpz_poly_coeff_division(mpz_poly * f, unsigned long m);

mpz_poly mpz_poly_coeff_product_si(mpz_poly * f, long m);

mpz_poly mpz_poly_sum_sorted(mpz_poly * f, mpz_poly * g);

mpz_poly mpz_poly_sum(mpz_poly * f, mpz_poly * g);

mpz_poly mpz_poly_prod(mpz_poly * f, mpz_poly * g);

int mpz_poly_is_constant(mpz_poly * f);

mpz_poly compose_psi(mpz_poly * f,long l,rootsieve_dictionary rdict);

mpz_poly compose_reduce( mpz_poly * f, long l,rootsieve_dictionary rdict,unsigned long m);

long rotation_inner(rootsieve_dictionary rdict, long p,mpz_poly * ff,mpz_poly * gg,mpz_poly * fdg_gdf,long u0,long v0,long l0,int ld,int m,long twist_v1,double scale,mpz_poly dphi);

void mpz_poly_eval_si(mpz_t r, mpz_poly * f, long a);

mpz_poly mpz_poly_identity();

void mpz_poly_constant_coeff(mpz_t op,mpz_poly * f);

long mpz_poly_constant_coeff_si(mpz_poly * f);

long modular_inverse( long r,  long N );

int fill_alpha(rootsieve_dictionary rdict, long p);

void extend_ppow_list(rootsieve_dictionary rdict,int l, long p);

void mpz_poly_set(mpz_poly * f,mpz_t * clist);

