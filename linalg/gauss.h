#define SAVE_KERNEL_MEMORY 1

extern int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
		  int limbs_per_row, int limbs_per_col);
