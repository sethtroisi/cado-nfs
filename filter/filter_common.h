#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#define STR(s) XSTR(s)
#define XSTR(s) #s

#define CA_DUP1 314159265358979323UL
#define CB_DUP1 271828182845904523UL

#define CA_DUP2 271828182845904523UL
#define CB_DUP2 577215664901532889UL

/* Maximum level for a merge. Such a large value is only useful when not using
 * BW. */
#define MERGE_LEVEL_MAX 256

#define DEFAULT_FILTER_EXCESS 160

#define DEFAULT_PURGE_NPASS 50
#define DEFAULT_PURGE_REQUIRED_EXCESS 0.1
#define DEFAULT_PURGE_NPT 4

#define DEFAULT_MERGE_MAXLEVEL 10
#define DEFAULT_MERGE_FORBW 0
#define DEFAULT_MERGE_RATIO 1.1
#define DEFAULT_MERGE_COVERNMAX 100.0
#define DEFAULT_MERGE_MKZTYPE 1 /* pure Markowitz */
#define DEFAULT_MERGE_WMSTMAX 7 /* relevant only if mkztype == 2 */
#ifndef FOR_DL
#define DEFAULT_MERGE_SKIP 32
#else
#define DEFAULT_MERGE_SKIP 0
#endif


#include "utils_with_io.h"

unsigned int insert_rel_in_table_no_e (earlyparsed_relation_ptr, index_t, index_t **, weight_t *);
unsigned int insert_rel_in_table_with_e (earlyparsed_relation_ptr, index_t,
                                         ideal_merge_t **, int32_t *);
int cmp_index (const void *p, const void *q);
int cmp_int (const void *p, const void *q);
int cmp_ideal_merge (const void *p, const void *q);
int cmp_int2 (const void *p, const void *q);

#endif /* FILTER_UTILS_H_ */
