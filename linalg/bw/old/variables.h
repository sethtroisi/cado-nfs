#ifndef VARIABLES_H_
#define VARIABLES_H_

#include "types.h"

#ifdef	__cplusplus
extern "C" {
#endif

extern coord_t nrows;
#ifdef HARDCODE_PARAMS
#define m_param HARD_m_param
#define n_param HARD_n_param
#define bigdim (m_param + n_param)
#else
#define m_param computed_m_param
#define n_param computed_n_param
#define bigdim computed_bigdim
#endif
extern int computed_m_param;
extern int computed_n_param;
extern int computed_bigdim;
extern int total_work;

extern void consistency_check(const char *, int, int);

#ifdef	__cplusplus
}
#endif

#endif /* VARIABLES_H_ */
