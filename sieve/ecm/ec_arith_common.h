#ifndef _EC_ARITH_COMMON_H_
#define _EC_ARITH_COMMON_H_

/* Types of coordinates */
typedef enum {
  SHORT_WEIERSTRASS_aff,
  SHORT_WEIERSTRASS_proj,
  MONTGOMERY_xz,
  TWISTED_EDWARDS_proj,
  TWISTED_EDWARDS_ext,
} ec_point_coord_type_t;


/* A point on an elliptic curve:
 *   The significant coordinates are:
 *      - x,y for SHORT_WEIERSTRASS_aff
 *      - x,y,z for SHORT_WEIERSTRASS_proj
 *      - x,z for MONTGOMERY_xz
 *      - x,y,z for TWISTED_EDWARDS_proj
 *      - x,y,z,t for TWISTED_EDWARDS_ext
 */
struct ec_point_s
{
  residue_t x,y,z,t;
};

typedef struct ec_point_s ec_point_t[1];
typedef struct ec_point_s *ec_point_ptr;
typedef const struct ec_point_s *ec_point_srcptr;

static inline void
ec_point_init (ec_point_t P, const modulus_t m)
{
  mod_init (P->x, m);
  mod_init (P->y, m);
  mod_init (P->z, m);
  mod_init (P->t, m);
}

static inline void
ec_point_clear (ec_point_t P, const modulus_t m)
{
  mod_clear (P->x, m);
  mod_clear (P->y, m);
  mod_clear (P->z, m);
  mod_clear (P->t, m);
}

static inline void
ec_point_set (ec_point_t Q, const ec_point_t P, const modulus_t m)
{
  mod_set (Q->x, P->x, m);
  mod_set (Q->y, P->y, m);
  mod_set (Q->z, P->z, m);
  mod_set (Q->t, P->t, m);
}

static inline void
ec_point_swap (ec_point_t Q, ec_point_t P, const modulus_t m)
{
  mod_swap (Q->x, P->x, m);
  mod_swap (Q->y, P->y, m);
  mod_swap (Q->t, P->t, m);
  mod_swap (Q->z, P->z, m);
}

static inline void
ec_point_fprintf (FILE *out, ec_point_t P, const ec_point_coord_type_t coord,
                  MAYBE_UNUSED const modulus_t m)
{
  // XXX experimental and maybe not optimal

  residue_t v1, v2, tmp;
  mod_init (v1, m);
  mod_init (v2, m);
  mod_init (tmp, m);
  
  mod_inv (v1, P->z, m);
  mod_mul (v1, v1, P->x, m);

  mod_sub (v2, P->z, P->y, m);
  mod_inv (v2, v2, m);
  mod_add (tmp, P->z, P->y, m);
  mod_mul (v2, v2, tmp, m);
  
#ifdef MOD_SIZE
  modint_t x, y, z, t, w1, w2;

  mod_intinit (x);
  mod_intinit (y);
  mod_intinit (z);
  mod_intinit (t);
  mod_intinit (w1);
  mod_intinit (w2);
  
  mod_get_int (x, P->x, m);
  mod_get_int (y, P->y, m);
  mod_get_int (z, P->z, m);
  mod_get_int (t, P->t, m);

  switch (coord)
  {
    case SHORT_WEIERSTRASS_aff:
      gmp_fprintf (out, "(%Nd, %Nd)", x, MOD_SIZE, y, MOD_SIZE);
      break;
    case SHORT_WEIERSTRASS_proj:
      gmp_fprintf (out, "(%Nd : %Nd : %Nd)", x, MOD_SIZE, y, MOD_SIZE, z,
                                                                      MOD_SIZE);
      break;
    case MONTGOMERY_xz:
      gmp_fprintf (out, "(%Nd :: %Nd)\n", x, MOD_SIZE, z, MOD_SIZE);
      mod_get_int (w1, v1, m);
      gmp_fprintf (out, "  X/Z: %Nd\n", w1, MOD_SIZE);
      break;
    case TWISTED_EDWARDS_proj:
      gmp_fprintf (out, "(%Nd : %Nd : %Nd)", x, MOD_SIZE, y, MOD_SIZE, z,
                                                                      MOD_SIZE);
      break;
    case TWISTED_EDWARDS_ext:
      gmp_fprintf (out, "(%Nd : %Nd : %Nd : %Nd)\n", x, MOD_SIZE, y, MOD_SIZE,
                                                   z, MOD_SIZE, t, MOD_SIZE);
      mod_get_int (w2, v2, m);
      gmp_fprintf (out, "  (Z+Y)/(Z-Y): %Nd\n", w2, MOD_SIZE);
      
      break;
    default:
      fprintf (stderr, "%s: unknown coordinates system\n", __func__);
      abort ();
  }

  mod_intclear (x);
  mod_intclear (y);
  mod_intclear (z);
  mod_intclear (t);
  mod_intclear (w1);
  mod_intclear (w2);
#else
  switch (coord)
  {
    case SHORT_WEIERSTRASS_aff:
      gmp_fprintf (out, "(%Zd, %Zd)", P->x, P->y);
      break;
    case SHORT_WEIERSTRASS_proj:
      gmp_fprintf (out, "(%Zd : %Zd : %Zd)", P->x, P->y, P->z);
      break;
    case MONTGOMERY_xz:
      gmp_fprintf (out, "(%Zd :: %Zd)\n", P->x, P->z);
      gmp_fprintf (out, "  X/Z = %Zd\n", v1);
      break;
    case TWISTED_EDWARDS_proj:
      gmp_fprintf (out, "(%Zd : %Zd : %Zd)", P->x, P->y, P->z);
      break;
    case TWISTED_EDWARDS_ext:
      gmp_fprintf (out, "(%Zd : %Zd : %Zd : %Zd)\n", P->x, P->y, P->z, P->t);
      gmp_fprintf (out, "  (Z+Y)/(Z-Y) = %Zd\n", v2);
      break;
    default:
      fprintf (stderr, "%s: unknown coordinates system\n", __func__);
      abort ();
  }
#endif
  mod_clear (v1, m);
  mod_clear (v2, m);
  mod_clear (tmp, m);
}

static inline void
ec_montgomery_curve_fprintf (FILE *out, residue_t A, ec_point_t P,
                             MAYBE_UNUSED const modulus_t m)
{
  // XXX experimental

#ifdef MOD_SIZE
  modint_t cc, mod;

  mod_intinit (cc);
  mod_intinit (mod);

  mod_get_int (cc, A, m);
  mod_getmod_int (mod, m);

  gmp_fprintf (out, "  Montgomery curve: B*Y^2 = X^3 + A*X^2*Z + X*Z^2\n"
	       "  A = %Nd\n  modulo %Nd\n", cc, MOD_SIZE, mod, MOD_SIZE);

  mod_intclear (cc);
  mod_intclear (mod);
#else
  gmp_fprintf (out, "  Montgomery curve: B*Y^2 = X^3 + A*X^2*Z + X*Z^2\n"
                    "  A = %Zd\n  modulo %Zd\n", A, m);

#endif

  if (P)
  {
    fputs ("  with point (X::Z) = ", out);
    ec_point_fprintf (out, P, MONTGOMERY_xz, m);
    fputc ('\n', out);
  }
}

static inline void
ec_twisted_edwards_ext_curve_fprintf (FILE *out, residue_t A, ec_point_t P,
                                      MAYBE_UNUSED const modulus_t m)
{
  // XXX experimental
#ifdef MOD_SIZE
  modint_t mod;
  mod_intinit (mod);
  mod_getmod_int (mod, m);
  gmp_fprintf (out, "Twisted Edwards curve: -X^2 + Y^2 = Z^2 + d*T^2\n"
	       "  XY = ZT (extended coordinates)\n  modulo %Nd\n", mod, MOD_SIZE);
  mod_intclear (mod);

#else
  gmp_fprintf (out, "Twisted Edwards curve: -X^2 + Y^2 = Z^2 + d*T^2\n"
                    "  XY = ZT (extended coordinates)\n  modulo %Zd\n", m);
#endif


  if (A)
  {
#ifdef MOD_SIZE
    modint_t cc;
    mod_intinit (cc);
    mod_get_int (cc, A, m);
    gmp_fprintf (out, "  Equivalent to Montgomery curve with\n"
		      "  A = %Nd\n", cc, MOD_SIZE);
    mod_intclear (cc);
#else
    gmp_fprintf (out, "  Equivalent to Montgomery curve with\n"
		      "  A = %Zd\n", A);
#endif
  }
  
  if (P)
    {
      fputs ("  with point (X:Y:Z:T) = ", out);
      ec_point_fprintf (out, P, TWISTED_EDWARDS_ext, m);
      fputc ('\n', out);
    }
  
}

static inline void
ec_weierstrass_curve_fprintf (FILE *out, residue_t a, ec_point_t P,
                             MAYBE_UNUSED const modulus_t m)
{
  // XXX experimental
#ifdef MOD_SIZE
  modint_t cc, mod;

  mod_intinit (cc);
  mod_intinit (mod);

  mod_get_int (cc, a, m);
  mod_getmod_int (mod, m);

  gmp_fprintf (out, "Weierstrass curve: y^2 = x^3 + a*x + b\n"
	       "  a = %Nd\n  modulo %Nd\n", cc, MOD_SIZE, mod, MOD_SIZE);

  mod_intclear (cc);
  mod_intclear (mod);
#else
  gmp_fprintf (out, "Weierstrass curve: y^2 = x^3 + a*x + b\n"
                    "  a = %Zd\n  modulo %Zd\n", a, m);
#endif

  if (P)
  {
    fputs ("  with point (x, y) = ", out);
    ec_point_fprintf (out, P, SHORT_WEIERSTRASS_aff, m);
    fputc ('\n', out);
  }
}

#endif /* _EC_ARITH_COMMON_H_ */
