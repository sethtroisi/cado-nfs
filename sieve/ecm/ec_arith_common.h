#ifndef _EC_ARITH_COMMON_H_
#define _EC_ARITH_COMMON_H_

/* Types of coordinates */
typedef enum {
  AFF,
  MONTG,
  EDW_proj,
  EDW_ext
} ec_point_coord_type_t;


/* A point on an elliptic curve */
/* In affine coordinates (AFF), only x and y are guaranteed */
/* In Montgomery coordinates (MONTG), only x and z are guaranteed */
/* In Edwards projective (EDW_proj), only x, y and z are guaranteed */
/* In Edwards extended (EDW_ext), x, y, z, t are guaranteed */
  
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
ec_point_print (ec_point_t P, const ec_point_coord_type_t coord_type)
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

#endif /* _EC_ARITH_COMMON_H_ */
