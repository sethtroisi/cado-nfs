import sys

HEAD = """\
#include "macros.h"

/* Typecasting inline functions that assert that casted values are within 
   range of the casted-to type.
   Format is cast_<source type>_<target type>()
   where each of <source type> and <target type> are one of
   ulong	for unsigned long int
   uint		for unsigned int
   long		for long int
   int		for int
   ushort	for unsigned short int
   short	for short int
   uchar	for unsigned char
   char		for char
   size		for size_t
   
   For each function, there is also a
   cast_<source type>_<target type>_fast()
   variant which asserts range validity only if WANT_ASSERT_CONVERSION is
   defined. Thus the non-fast variant should be used in code that is not
   speed critical, and the fast variant in code that is.
*/

#ifdef WANT_ASSERT_CONVERSION
#define ASSERT_CONVERSION(x) assert(x)
#else
#define ASSERT_CONVERSION(x)
#endif

"""

SRC_NAMES = ("srctype", "srctype_c", "srcmin", "srcmax")
DST_NAMES = ("dsttype", "dsttype_c", "dstmin", "dstmax")
TYPES = (
  ("ulong", "unsigned long int", "0", "ULONG_MAX"),
  ("uint", "unsigned int", "0", "UINT_MAX"),
  ("long", "long int", "LONG_MIN", "LONG_MAX"),
  ("int", "int", "INT_MIN", "INT_MAX"),
  ("ushort", "unsigned short int", "0", "USHRT_MAX"), # The 'O' would have made it excessively long?
  ("short", "short int", "SHRT_MIN", "SHRT_MAX"),
  ("uchar", "unsigned char", "0", "UCHAR_MAX"), # The 'O' would have made it excessively long?
  ("char", "char", "CHAR_MIN", "CHAR_MAX"),
  ("size", "size_t", "0", "SIZE_MAX")
)

TEMPLATE1 = """\
static inline {dsttype_c}
cast_{srctype}_{dsttype}{fast}({srctype_c} x)
{{
"""
TEMPLATE2 = """\
#if {srcmin} < {dstmin}
    {assert}(x >= {dstmin});
#endif
"""
TEMPLATE3 = """\
#if {srcmax} > {dstmax}
    {assert}(x <= {dstmax});
#endif
"""
TEMPLATE4 = """\
  return ({dsttype_c}) x;
}}
"""

SLOW = [("assert", "ASSERT"), ("fast", "")]
FAST = [("assert", "ASSERT_CONVERSION"), ("fast", "_fast")]

def print_template(d):
  sys.stdout.write(TEMPLATE1.format(**d))
  if not d["srcmin"] == d["dstmin"]:
    sys.stdout.write(TEMPLATE2.format(**d))
  if not d["srcmax"] == d["dstmax"]:
    sys.stdout.write(TEMPLATE3.format(**d))
  sys.stdout.write(TEMPLATE4.format(**d))
  sys.stdout.write("\n")

print(HEAD)
for srctype in TYPES:
  for dsttype in TYPES:
    if srctype is dsttype:
      continue
    values = list(zip(SRC_NAMES, srctype)) + list(zip(DST_NAMES, dsttype))
    print_template(dict(values + SLOW))
    print_template(dict(values + FAST))
