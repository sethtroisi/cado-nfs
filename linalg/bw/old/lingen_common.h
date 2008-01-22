#ifndef LINGEN_COMMON_H_
#define LINGEN_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "structure.h"

struct col_id {
	int degnom;
	int pos;
};

struct bw_context {
	bw_mnpoly	a;
	int		t0;
	int		cur_deg;	/* actually the number of coeffs
					   available ; that is, the degree
					   reached + 1 */
};

struct bw_iterator {
	bw_nbpoly	f;
	bw_mbmat	ctaf;
	int		t;
	struct col_id * clist;
	int	      * pivots_list;
	unsigned int  * chance_list;
};

#define	Lmacro(N, m, n)	(iceildiv((N)+2*(n),(m))+iceildiv((N)+2*(n),(n))+10)

typedef int (*sortfunc_t)(const void*, const void*);

extern void bw_commit_f(bw_nbpoly fpoly, const int *);
extern void print_chance_list(unsigned int t, unsigned int *);
extern int read_data_for_series(bw_mnpoly a);
extern void give_mnpoly_rank_info(bw_mnpoly a, int deg);
extern void read_mat_file_header(const char * name);

extern const char * a_meta_filename;
extern const char * valu_meta_filename;
extern const char * pi_meta_filename;

#ifdef __cplusplus
}
#endif

#endif	/* LINGEN_COMMON_H_ */
