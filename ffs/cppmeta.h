#ifndef __CPPMETA_H__
#define __CPPMETA_H__



// Concatenating macro, which first expands its arguments.
#define __CAT(x, y) x##y
#define   CAT(x, y) __CAT(x, y)


// Expand a tuple.
#define __EXPAND(...) __VA_ARGS__
#define   EXPAND(l)   __EXPAND l


// Count the number of arguments (up to 4).
#define __NUM_ARGS(x4, x3, x2, x1, n, ...) n
#define   NUM_ARGS(...) __NUM_ARGS(__VA_ARGS__, 4, 3, 2, 1, 0, )


// Apply macro <f> to a list of arguments.
#define __FOR_ALL_0(f, a)
#define __FOR_ALL_1(f, a, x)      f(a, x, 1)
#define __FOR_ALL_2(f, a, x, ...) f(a, x, 2), __FOR_ALL_1(f, a, __VA_ARGS__)
#define __FOR_ALL_3(f, a, x, ...) f(a, x, 3), __FOR_ALL_2(f, a, __VA_ARGS__)
#define __FOR_ALL_4(f, a, x, ...) f(a, x, 4), __FOR_ALL_3(f, a, __VA_ARGS__)
#define __FOR_ALL2(n, f, a, ...)  __FOR_ALL_##n(f, a, __VA_ARGS__)
#define __FOR_ALL1(n, f, a, ...)  __FOR_ALL2(n, f, a, __VA_ARGS__)
#define   FOR_ALL(f, a, ...)      __FOR_ALL1(NUM_ARGS(__VA_ARGS__), \
                                             f, a, __VA_ARGS__)


// Apply macro <f> to integers n-1 downto 0 (n up to 16).
#define __FOR_0( f, a)
#define __FOR_1( f, a) f(a,  0)
#define __FOR_2( f, a) f(a,  1) __FOR_1 (f, a)
#define __FOR_3( f, a) f(a,  2) __FOR_2 (f, a)
#define __FOR_4( f, a) f(a,  3) __FOR_3 (f, a)
#define __FOR_5( f, a) f(a,  4) __FOR_4 (f, a)
#define __FOR_6( f, a) f(a,  5) __FOR_5 (f, a)
#define __FOR_7( f, a) f(a,  6) __FOR_6 (f, a)
#define __FOR_8( f, a) f(a,  7) __FOR_7 (f, a)
#define __FOR_9( f, a) f(a,  8) __FOR_8 (f, a)
#define __FOR_10(f, a) f(a,  9) __FOR_9 (f, a)
#define __FOR_11(f, a) f(a, 10) __FOR_10(f, a)
#define __FOR_12(f, a) f(a, 11) __FOR_11(f, a)
#define __FOR_13(f, a) f(a, 12) __FOR_12(f, a)
#define __FOR_14(f, a) f(a, 13) __FOR_13(f, a)
#define __FOR_15(f, a) f(a, 14) __FOR_14(f, a)
#define __FOR_16(f, a) f(a, 15) __FOR_15(f, a)
#define __FOR(n, f, a) __FOR_##n(f, a)
#define   FOR(f, a, n) __FOR(n, f, a)


// Return <then> if __<var>_<cond> is defined to "," (a single comma),
// otherwise return <else>.
#define __IF2(x1, x2, ret, ...)     ret
#define __IF1(...)                  __IF2(__VA_ARGS__)
#define   IF(var, cond, then, else) __IF1(CAT(CAT(__, var), CAT(_, cond)), \
                                          then, else, )

// Return <val> if __<var>_<cond> is defined to ", <val>" (a comma followed
// by <val>), otherwise return <default>.
#define SWITCH(var, cond, default)  __IF1(, CAT(CAT(__, var), CAT(_, cond)), \
                                          default, )


// Useful conditions.
#undef  EMPTY
#define ___EMPTY       ,

#undef  EQ
#undef  LE
#define __CMP_2_2_EQ   ,
#define __CMP_16_16_EQ ,
#define __CMP_32_32_EQ ,
#define __CMP_64_64_EQ ,
#define __CMP_16_32_LE ,
#define __CMP_16_64_LE ,
#define __CMP_32_64_LE ,
#define   CMP(x, y) CAT(CAT(CMP_, x), CAT(_, y))

#endif  /* __CPPMETA_H__ */
