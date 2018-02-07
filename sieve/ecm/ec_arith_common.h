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
                  const modulus_t m)
{
  modint_t x, y, z, t;

  mod_intinit (x);
  mod_intinit (y);
  mod_intinit (z);
  mod_intinit (t);
  
  mod_get_int (x, P->x, m);
  mod_get_int (y, P->y, m);
  mod_get_int (z, P->z, m);
  mod_get_int (t, P->t, m);

  switch (coord)
  {
    case SHORT_WEIERSTRASS_aff: /* (x, y) */
      mod_fprintf (out, "%s(0x%" PRIMODx ", 0x%" PRIMODx ")", MOD_PRINT_INT(x),
          MOD_PRINT_INT(y));
      break;
    case TWISTED_EDWARDS_proj: /* (X : Y : Z) */
    case SHORT_WEIERSTRASS_proj: /* (X : Y : Z) */
      mod_fprintf (out, "(0x%" PRIMODx " : 0x%" PRIMODx " : 0x%" PRIMODx ")",
          MOD_PRINT_INT(x), MOD_PRINT_INT(y), MOD_PRINT_INT(z));
      break;
    case MONTGOMERY_xz: /* (X :: Z) */
      mod_fprintf (out, "(0x%" PRIMODx " :: 0x%" PRIMODx ")", MOD_PRINT_INT(x),
          MOD_PRINT_INT(z));
      break;
    case TWISTED_EDWARDS_ext: /* (X : Y : Z : T) */
      mod_fprintf (out, "(0x%" PRIMODx " : 0x%" PRIMODx " : 0x%" PRIMODx " : "
          "0x%" PRIMODx ")", MOD_PRINT_INT(x), MOD_PRINT_INT(y),
          MOD_PRINT_INT(z), MOD_PRINT_INT(t));
      break;
    default:
      fprintf (stderr, "%s: unknown coordinates system\n", __func__);
      abort ();
  }

  mod_intclear (x);
  mod_intclear (y);
  mod_intclear (z);
  mod_intclear (t);
}

static inline void
ec_montgomery_curve_fprintf (FILE *out, const char *prefix, residue_t A,
                             ec_point_t P, const modulus_t m)
{
  modint_t cc;
  const char *pre = (prefix == NULL) ? "" : prefix;

  mod_intinit (cc);
  mod_get_int (cc, A, m);

  mod_fprintf (out, "%sMontgomery curve: B*Y^2 = X^3 + A*X^2*Z + X*Z^2\n"
                    "%sA = 0x%" PRIMODx "\n", pre, pre, MOD_PRINT_INT (cc));

  mod_intclear (cc);

  if (P)
  {
    fprintf (out, "%swith point (X::Z) = ", pre);
    ec_point_fprintf (out, P, MONTGOMERY_xz, m);
    fputc ('\n', out);
  }
}

/* A is the curve coefficient of the isomorphic Montgomery curve, the relation
 * is (because we have a=-1):
 *    A = 2*(1-d)/(1+d)
 *    d = -(A-2)/(A+2)
 */
static inline void
ec_twisted_edwards_ext_curve_fprintf (FILE *out, const char *prefix,
                                      residue_t d, ec_point_t P,
                                      const modulus_t m)
{
  const char *pre = (prefix == NULL) ? "" : prefix;

  modint_t cc;
  mod_intinit (cc);

  mod_get_int (cc, d, m);

  mod_fprintf (out, "%sTwisted Edwards curve: -X^2 + Y^2 = Z^2 + d*T^2\n"
                    "%sXY = ZT (extended coordinates)\n%sd = 0x%" PRIMODx "\n",
                    pre, pre, pre, MOD_PRINT_INT (cc));

  mod_intclear (cc);

  if (P)
  {
    fprintf (out, "%swith point (X:Y:Z:T) = ", pre);
    ec_point_fprintf (out, P, TWISTED_EDWARDS_ext, m);
    fputc ('\n', out);
  }
}

static inline void
ec_weierstrass_curve_fprintf (FILE *out, residue_t a, ec_point_t P,
                             MAYBE_UNUSED const modulus_t m)
{
  // XXX experimental
  // TODO rewrite it as the others
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
