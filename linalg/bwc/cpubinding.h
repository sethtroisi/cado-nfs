#ifndef CPUBINDING_H_
#define CPUBINDING_H_

#ifdef __cplusplus
extern "C" {
#endif

/* This returns an opaque pointer to data which will be used to perform
 * the actual cpu binding. This function must be called in
 * single-threaded context.
 * 
 * This returns NULL if cpubinding failed.
 */
void * cpubinding_get_info(char ** messages, param_list_ptr pl, int thread_split[2]);

/* perform the actual pinning. This must be called for each thread */
void cpubinding_do_pinning(void * pinning_info, int i, int j);

/* free the opaque pointer */
void cpubinding_free_info(void * pinning_info, int thread_split[2]);

#ifdef __cplusplus
}
#endif

#endif	/* CPUBINDING_H_ */
