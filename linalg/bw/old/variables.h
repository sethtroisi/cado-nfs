#ifndef VARIABLES_H_
#define VARIABLES_H_

#ifdef	__cplusplus
extern "C" {
#endif

extern coord_t nrows;
extern coord_t ncols;
#ifdef HARDCODE_PARAMS
#define m_param HARD_m_param
#define n_param HARD_n_param
#else
#define m_param computed_m_param
#define n_param computed_n_param
#endif
extern int computed_m_param;
extern int computed_n_param;
extern int bank_num;
extern int total_work;
extern int periodicity;

#define subdim	n_param
#define bigdim	(subdim + m_param)

extern void consistency_check(const char *, int, int);

#ifdef	__cplusplus
}
#endif

#endif /* VARIABLES_H_ */
