#ifndef SINGLETON_REMOVAL_H_
#define SINGLETON_REMOVAL_H_

void singleton_removal_oneiter_mono (purge_matrix_ptr mat);
void singleton_removal_oneiter_mt (purge_matrix_ptr, unsigned int);
int64_t singleton_removal (purge_matrix_ptr, unsigned int, int);


#endif /* SINGLETON_REMOVAL_H_ */
