#ifndef _ELL_ARITH_H_
#define _ELL_ARITH_H_

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

void ell_point_init (ell_point_t, const modulus_t);
void ell_point_clear (ell_point_t, const modulus_t);
void ell_point_set (ell_point_t, const ell_point_t, const modulus_t);
void ell_point_swap (ell_point_t, ell_point_t, const modulus_t);
void ell_point_print (ell_point_t, const ell_point_coord_type_t);

#endif
