#ifndef __GRAY_H__
#define __GRAY_H__

#include "cppmeta.h"



// Binary Gray code.
#if defined(USE_F2)
# define __GRAY_1            0
# define __GRAY_2  __GRAY_1, 1, __GRAY_1
# define __GRAY_3  __GRAY_2, 2, __GRAY_2
# define __GRAY_4  __GRAY_3, 3, __GRAY_3
# define __GRAY_5  __GRAY_4, 4, __GRAY_4
# define __GRAY_6  __GRAY_5, 5, __GRAY_5
# define __GRAY_7  __GRAY_6, 6, __GRAY_6
# define __GRAY_8  __GRAY_7, 7, __GRAY_7
# define __GRAY_9  __GRAY_8, 8, __GRAY_8
# define __GRAY_10 __GRAY_9, 9, __GRAY_9

  // Length of size-n Gray code is 2^n-1.
# define GRAY_LENGTH(n) ((1u<<(n))-1)

// Ternary Gray code.
#elif defined(USE_F3)
# define __GRAY_1            0,           0
# define __GRAY_2  __GRAY_1, 1, __GRAY_1, 1, __GRAY_1
# define __GRAY_3  __GRAY_2, 2, __GRAY_2, 2, __GRAY_2
# define __GRAY_4  __GRAY_3, 3, __GRAY_3, 3, __GRAY_3
# define __GRAY_5  __GRAY_4, 4, __GRAY_4, 4, __GRAY_4
# define __GRAY_6  __GRAY_5, 5, __GRAY_5, 5, __GRAY_5
# define __GRAY_7  __GRAY_6, 6, __GRAY_6, 6, __GRAY_6
# define __GRAY_8  __GRAY_7, 7, __GRAY_7, 7, __GRAY_7
# define __GRAY_9  __GRAY_8, 8, __GRAY_8, 8, __GRAY_8
# define __GRAY_10 __GRAY_9, 9, __GRAY_9, 9, __GRAY_9

  // Length of size-n Gray code is 3^n-1.
  static const unsigned __f3_gray_length[] = {
    2, 8, 26, 80, 242, 728, 2186, 6560, 19682, 59048
  };
# define GRAY_LENGTH(n) __f3_gray_length[(n)-1]

#endif

#define GRAY(n) CAT(__GRAY_, n)

#endif  /* __GRAY_H__ */
