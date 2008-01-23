#ifndef SLAVE_H_
#define SLAVE_H_

#ifdef	__cplusplus
extern "C" {
#endif

extern int total_work;
extern int periodicity;

extern coord_t *bw_x;
extern stype32 * bw_m;
extern stype32 * bw_x0;

#ifdef HARDCODE_PARAMS
#define nbxs HARD_nbxs
#define nbys HARD_nbys
#else
extern int computed_nbxs;
extern int computed_nbys;
#define nbxs computed_nbxs
#define nbys computed_nbys
#endif

#ifdef CORNERS_H_
extern corner_s first,last;
#endif

enum task_type_t { TASK_DOTPRODUCTS, TASK_POLYNOMIAL };

struct mksol_info_block {
	int	solution_col;
	int	t_value;
	int	valuation;
	void  * sum;
	int	differential;
};

extern enum task_type_t task_type;

#ifdef	__cplusplus
}
#endif

#endif /* SLAVE_H_ */
