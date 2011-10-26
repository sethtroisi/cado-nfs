#ifndef LAS_UNSIEVE_H_
#define LAS_UNSIEVE_H_

struct unsieve_aux_data_s {
    unsigned int *lpf; /* lpf[i] is largest prime factor of i, for i < I */
};
typedef struct unsieve_aux_data_s unsieve_aux_data[1];
typedef struct unsieve_aux_data_s * unsieve_aux_data_ptr;

#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

void sieve_info_init_unsieve_data(sieve_info_ptr si);
void sieve_info_clear_unsieve_data(sieve_info_ptr si);
void unsieve_not_coprime (unsigned char *S, const int N, sieve_info_srcptr si);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_UNSIEVE_H_ */
