#ifndef BYTECODE_MISHMASH_B1_DATA_H
#define BYTECODE_MISHMASH_B1_DATA_H

struct mishmash_B1_data_s
{
  unsigned int B1;
  unsigned int len; /* length of the bytecode */
  bytecode_const bc;
};

typedef struct mishmash_B1_data_s mishmash_B1_data_t;

static const uint8_t _B1_5_bc[] = {
    0x00,
    0x11, 0x71, 0x04, /* 3*5 */
    0xff,
  }; /* end of bytecode for B1=5 */
static const uint8_t _B1_105_bc[] = {
    0x03,
    0x11, 0xf1, 0x0c, 0x0c, /* 97*43*37*31*13*7*5 */
    0x11, 0xf1, 0x08, 0x0c, /* 73*71*61*17*5 */
    0x11, 0xc1, 0x01, 0x0b, 0x71, 0x09, /* 89*53*29*23 */
    0x21, 0xa2, 0x01, 0xff, 0x12, 0x51, 0x11, 0x62, 0x05, /* 101*83*79*19 */
    0x21, 0xa2, 0x05, 0x80, 0x02, 0xff, 0x10, 0x41, 0x07, 0x72, 0x04, /* 103*67*59*11 */
    0x81, 0x69, 0x66, 0x69, 0x66, 0x69, 0x66, 0x69, 0x66, 
          0x69, 0X04, 0x66, 
          0x69, 0x0b, 0x0b, 0x03, 0x03, 0x0b, 0x03, 0x66,
          0x69, 0x0b, 0x0b, 0x0b, 0x0b, 0x03, 0x03, 0x46, /* 3^4*7*41*47 */
    0xff,
  }; /* end of bytecode for B1=105 */

#define mishmash_B1_data_len (sizeof(mishmash_B1_data)/sizeof(mishmash_B1_data_t))
/* Sorted by increasing B1 */
static mishmash_B1_data_t mishmash_B1_data[] = {
    (mishmash_B1_data_t) { .B1 = 5,
                           .len = sizeof(_B1_5_bc)/sizeof(*_B1_5_bc),
                           .bc = _B1_5_bc },
    (mishmash_B1_data_t) { .B1 = 105,
                           .len = sizeof(_B1_105_bc)/sizeof(*_B1_105_bc),
                           .bc = _B1_105_bc },
  }; /* end of B1_data array */

#endif /* BYTECODE_MISHMASH_B1_DATA_H */

