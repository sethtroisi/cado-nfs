/* Data type to store column values: a 32-bit integer should be enough for
   most applications!
*/
#include <stdint.h>
#define INT int32_t

extern void fprintRow(FILE *file, INT *row);
extern void removeWeight(INT **rows, int *wt, int i);
extern void addWeight(INT **rows, int *wt, int i);
extern void addRows(INT **rows, int i1, int i2, int len0);
