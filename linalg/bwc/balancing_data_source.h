#ifndef DATA_SOURCE_H_
#define DATA_SOURCE_H_

#include <stddef.h>
#include <stdint.h>

struct data_source_s {/*{{{*/
    size_t (*get)(void*,uint32_t **,size_t);
    size_t pos;
};
typedef struct data_source_s data_source[1];
typedef struct data_source_s *data_source_ptr;

/*}}}*/

#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif

#endif	/* DATA_SOURCE_H_ */
