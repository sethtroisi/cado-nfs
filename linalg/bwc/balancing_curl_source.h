#ifndef CURL_SOURCE_H_
#define CURL_SOURCE_H_

#include <curl/curl.h>
/* #include <curl/types.h> */
#include <curl/easy.h>

#include "balancing_data_source.h"
#include "ringbuf.h"

struct curl_source_s {
    data_source b;
    ringbuf r;
    CURL *curl_handle;
    const char * uri;
    pthread_t producer;
};

typedef struct curl_source_s curl_source[1];
typedef struct curl_source_s * curl_source_ptr;
typedef const struct curl_source_s * curl_source_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern data_source_ptr curl_source_alloc(const char * uri, size_t esz);
extern void curl_source_free(data_source_ptr q);
extern size_t curl_source_get(curl_source_ptr f, void ** p, size_t avail);

#ifdef __cplusplus
}
#endif

#endif	/* CURL_SOURCE_H_ */
