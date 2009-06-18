#ifndef CRC_H_
#define CRC_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This computes the crc32 value of the n 32-bit values pointed to by c.
 * XXX Having c of type unsigned long * is clearly a bug.
 */
extern uint32_t crc32(unsigned long * c, int n);

#ifdef __cplusplus
}
#endif

#endif	/* CRC_H_ */
