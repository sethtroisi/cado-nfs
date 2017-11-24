#ifndef _ECM_ARITH_COMMON_H_
#define _ECM_ARITH_COMMON_H_

/* Types of coordinates */
typedef enum {
  AFF,
  MONTG,
  EDW_proj,
  EDW_ext
} ell_point_coord_type_t;


/* A point on an elliptic curve */
/* In affine coordinates (AFF), only x and y are guaranteed */
/* In Montgomery coordinates (MONTG), only x and z are guaranteed */
/* In Edwards projective (EDW_proj), only x, y and z are guaranteed */
/* In Edwards extended (EDW_ext), x, y, z, t are guaranteed */
  
struct ell_point_s
{
  residue_t x,y,z,t;
};

typedef struct ell_point_s ell_point_t[1];
typedef struct ell_point_s *ell_point_ptr;
typedef const struct ell_point_s *ell_point_srcptr;

static inline void
ell_point_init (ell_point_t P, const modulus_t m)
{
  mod_init (P->x, m);
  mod_init (P->y, m);
  mod_init (P->z, m);
  mod_init (P->t, m);
}

static inline void
ell_point_clear (ell_point_t P, const modulus_t m)
{
  mod_clear (P->x, m);
  mod_clear (P->y, m);
  mod_clear (P->z, m);
  mod_clear (P->t, m);
}

static inline void
ell_point_set (ell_point_t Q, const ell_point_t P, const modulus_t m)
{
  mod_set (Q->x, P->x, m);
  mod_set (Q->y, P->y, m);
  mod_set (Q->z, P->z, m);
  mod_set (Q->t, P->t, m);
}

static inline void
ell_point_swap (ell_point_t Q, ell_point_t P, const modulus_t m)
{
  mod_swap (Q->x, P->x, m);
  mod_swap (Q->y, P->y, m);
  mod_swap (Q->t, P->t, m);
  mod_swap (Q->z, P->z, m);
}

static inline void
ell_point_print (ell_point_t P, const ell_point_coord_type_t coord_type)
  {
  /* FIXME need multiple precision print */
  
  printf ("(%lu", mod_intget_ul(P->x));
  if (coord_type == MONTG)
    printf (": %lu)", mod_intget_ul(P->z));
  else {
    printf (": %lu", mod_intget_ul(P->y));
    switch (coord_type)
      {
      case EDW_proj :
	printf (" : %lu", mod_intget_ul(P->z));
	break;
      case EDW_ext :
	printf (" : %lu : %lu", mod_intget_ul(P->t), mod_intget_ul(P->z));
	break;
      default :
	printf (")\n");
      }
  }
}

#endif
