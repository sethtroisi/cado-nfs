

struct ell_point_s
{
  residue_t x,y,z,t;
};

typedef struct ell_point_s ell_point_t[1];
typedef struct ell_point_s *ell_point_ptr;
typedef const struct ell_point_s *ell_point_srcptr;

enum ell_point_coord_type {AFF, MONTG, EDW_proj, EDW_ext} ;
